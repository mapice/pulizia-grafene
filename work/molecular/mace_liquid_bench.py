"""Bounded CPU/Metal feasibility test of the qualified molecular model.

Removing additive, element-specific atomic reference energies is an exact
Hamiltonian gauge for a fixed-composition system; forces/stress unchanged.
This prevents loss of small liquid-energy changes in float32 cancellation.
No density/dielectric/cleaning qualification follows from a speed benchmark.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,hashlib,time,traceback,gc
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
from ase import Atoms
from ase.io import read
from mace.calculators import MACECalculator
HERE=Path(__file__).resolve().parent
SHA='fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'


def structure(name):
    if name=='monomer':
        a=read(HERE/'results/clusters/HFIP.refined.xyz')
        a.positions-=a.positions.mean(0);return a
    p=HERE/'results/gaff-systems/HFIP_n0_seed20261008_N64'
    status=json.loads((p/'run-status.json').read_text())
    assert status[-1]['stage']=='production' and status[-1]['exitcode']==0
    # ASE reads native GRO element labels imperfectly; atom order is
    # independently tied to the original twelve-atom HFIP molecule.
    mono=read(HERE/'results/clusters/HFIP.refined.xyz')
    gro=read(p/'production.gro')
    count=8 if name=='liquid8' else 64
    a=Atoms(numbers=np.tile(mono.numbers,count),positions=gro.positions[:count*12],cell=gro.cell,pbc=True)
    assert len(a)==12*count and min(a.cell.lengths())>20
    return a


def main():
    p=argparse.ArgumentParser();p.add_argument('--device',choices=['cpu','mps'],default='cpu')
    p.add_argument('--dtype',choices=['float64','float32'],default='float64')
    p.add_argument('--structure',choices=['monomer','liquid8','liquid64'],default='monomer')
    p.add_argument('--tag',default='official');p.add_argument('--edge-chunk',type=int,default=0)
    p.add_argument('--fused-chunk',type=int,default=0)
    p.add_argument('--node-chunk',type=int,default=0)
    p.add_argument('--recompute-interactions',action='store_true')
    p.add_argument('--radial-chunk',type=int,default=0);args=p.parse_args()
    assert not (args.edge_chunk and args.fused_chunk)
    assert not (args.radial_chunk and (args.edge_chunk or args.fused_chunk or args.recompute_interactions or args.tag.startswith('checkpoint')))
    # Different benchmark names can select the explicit memory experiment.
    checkpointed=args.tag.startswith('checkpoint')
    assert not(args.device=='mps' and args.dtype=='float64')
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.set_default_dtype(torch.float64 if args.dtype=='float64' else torch.float32)
    if args.device=='mps':
        # Bound this process's unified GPU allocation; no OS/global setting.
        torch.mps.set_per_process_memory_fraction(.25)
    torch.serialization.add_safe_globals([slice])
    model=torch.load(path,map_location='cpu',weights_only=False)
    model=model.double() if args.dtype=='float64' else model.float()
    wrapped=[]
    if checkpointed:
        from memory_checkpoint import wrap_model
        wrapped=wrap_model(model)
    if args.edge_chunk:
        from chunk_tensor_product import wrap_products
        wrapped+=wrap_products(model,args.edge_chunk)
    if args.fused_chunk:
        from fused_chunk_message import wrap_messages
        wrapped+=wrap_messages(model,args.fused_chunk)
    if args.node_chunk:
        from chunk_node_products import wrap_node_products
        wrapped+=wrap_node_products(model,args.node_chunk)
    if args.recompute_interactions:
        from recompute_interaction import wrap_interactions
        wrapped+=wrap_interactions(model)
    if args.radial_chunk:
        from fused_radial_interaction import wrap_radial_interactions
        wrapped+=wrap_radial_interactions(model,args.radial_chunk)
    references=model.atomic_energies_fn.atomic_energies.detach().cpu().numpy().copy()
    z=[int(x) for x in model.atomic_numbers]
    a=structure(args.structure)
    a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    atomic_reference=sum(float(references.reshape(-1,len(z))[0,z.index(int(n))]) for n in a.numbers)
    root=HERE/'results/neural-liquid-bench';root.mkdir(exist_ok=True)
    label=f'{args.structure}-{args.device}-{args.dtype}-{args.tag}'
    target=root/(label+'.json')
    if target.exists():raise RuntimeError('Existing benchmark preserved')
    result=dict(model='MACE-POLAR-1-M',model_sha256=SHA,device=args.device,dtype=args.dtype,
                atoms=len(a),periodic=bool(a.pbc.all()),cell_A=a.cell.array.tolist(),
                input_sha256=hashlib.sha256(a.positions.tobytes()+a.numbers.tobytes()+a.cell.array.tobytes()).hexdigest(),
                torch_version=torch.__version__,scope='Feasibility/numerical benchmark, not liquid or cleaning qualification')
    result['installed_module_sha256']={name:hashlib.sha256((HERE/'.venv/lib/python3.12/site-packages/mace/modules'/name).read_bytes()).hexdigest()
                                     for name in ['extensions.py','field_blocks.py']}
    result['checkpointed_modules']=wrapped
    result['edge_chunk_size']=args.edge_chunk
    result['fused_chunk_size']=args.fused_chunk
    result['node_chunk_size']=args.node_chunk
    result['recompute_interactions']=args.recompute_interactions
    result['radial_chunk_size']=args.radial_chunk
    result['structure_note']='Dilute periodic numerical control, not a liquid' if args.structure=='liquid8' else args.structure
    try:
        # Check the gauge on CPU64 before changing any atomic references.
        if args.device=='cpu' and args.dtype=='float64' and args.structure=='monomer':
            c=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')
            b=a.copy();b.calc=c;original=float(b.get_potential_energy());original_force=b.get_forces().copy()
            with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
            c.reset();b=a.copy();b.calc=c;shifted=float(b.get_potential_energy());shifted_force=b.get_forces()
            error=original-shifted-atomic_reference
            force_error=float(abs(original_force-shifted_force).max())
            assert abs(error)<1e-8 and force_error<1e-10,(error,force_error)
            result['exact_gauge_check']=dict(offset_eV=atomic_reference,energy_identity_error_eV=error,
                                             force_difference_eV_A=force_error)
            del c,b;gc.collect()
        with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
        calc=MACECalculator(models=model,model_type='PolarMACE',device=args.device,default_dtype=args.dtype)
        a.calc=calc
        started=time.monotonic();energy=float(a.get_potential_energy());force=a.get_forces()
        stress=a.get_stress() if a.pbc.all() else None
        if args.device=='mps':torch.mps.synchronize()
        first=time.monotonic()-started
        times=[]
        for i in range(3):
            b=a.copy();b.positions[0,0]+=1e-4*(i+1);b.calc=calc
            started=time.monotonic();b.get_potential_energy();b.get_forces()
            if args.device=='mps':torch.mps.synchronize()
            times.append(time.monotonic()-started)
        charges=np.asarray(calc.results['charges'])
        result.update(accepted_finite_evaluation=True,shifted_energy_eV=energy,
                      atomic_reference_offset_eV=atomic_reference,forces_eV_A=force.tolist(),
                      stress_eV_A3=None if stress is None else np.asarray(stress).tolist(),
                      first_evaluation_seconds=first,median_force_seconds=float(np.median(times)),
                      net_charge_e=float(charges.sum()),
                      step1fs_projected_ps_per_hour=3.6/float(np.median(times)),
                      projection_is_not_an_executed_trajectory=True)
        if args.device=='mps':
            result['MPS_current_allocated_bytes']=torch.mps.current_allocated_memory()
            result['MPS_driver_allocated_bytes']=torch.mps.driver_allocated_memory()
            result['MPS_memory_fraction_limit']=.25
        assert np.isfinite(energy) and np.isfinite(force).all()
    except Exception as exc:
        result.update(accepted_finite_evaluation=False,error=str(exc),traceback=traceback.format_exc())
        if args.device=='mps':
            result.update(MPS_current_allocated_bytes=torch.mps.current_allocated_memory(),
                          MPS_driver_allocated_bytes=torch.mps.driver_allocated_memory(),
                          MPS_memory_fraction_limit=.25)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['forces_eV_A','traceback']},flush=True)
    if not result['accepted_finite_evaluation']:raise SystemExit(1)


if __name__=='__main__':main()
