"""Prepare small periodic cells for real, affordable liquid-model checks.

Gas monomer relaxed with the unchanged model, then periodic PackMol.
Experimental density is an INITIAL CONDITION, never a calculated result.
Several cell sizes and adequate sampling remain mandatory.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,hashlib,subprocess
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.set_default_dtype(torch.float64)
from ase.io import read,write
from ase.optimize import FIRE
from mace.calculators import MACECalculator
from neural_structure import SHA

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--molecules',type=int,choices=[16,32,64],default=16)
    p.add_argument('--seed',type=int,default=20261016)
    p.add_argument('--density',type=float,default=1.607);args=p.parse_args()
    assert .8<=args.density<=1.8
    suffix='' if abs(args.density-1.607)<1e-9 else f'_rho{args.density:.3f}'
    root=HERE/'results/neural-bulk'/(f'HFIP_N{args.molecules}_seed{args.seed}'+suffix)
    root.mkdir(parents=True,exist_ok=True)
    assert not (root/'prepared.json').exists()
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).double()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')
    a=read(HERE/'results/clusters/HFIP.refined.xyz');a.positions-=a.positions.mean(0)
    a.info.update(charge=0,spin=1,external_field=[0.,0.,0.]);a.calc=calc
    opt=FIRE(a,logfile=str(root/'gas-minimum.txt'),maxstep=.08)
    converged=bool(opt.run(fmax=.015,steps=200));assert converged
    fmax=float(np.linalg.norm(a.get_forces(),axis=1).max())
    # XYZ atom order is fixed and carries no implicit topology assumption.
    monomer=a.copy();monomer.positions-=monomer.positions.mean(0);monomer.calc=None
    write(root/'monomer.xyz',monomer,format='xyz')
    molar_mass=float(a.get_masses().sum());rho=args.density
    length=(args.molecules*molar_mass/(.602214076*rho))**(1/3)
    rmax=float(model.r_max);assert length>2*rmax+.2
    inp=f'''tolerance 2.0
filetype xyz
output packed.xyz
seed {args.seed}
nloop 2000
pbc {length:.10f} {length:.10f} {length:.10f}
structure monomer.xyz
  number {args.molecules}
end structure
'''
    (root/'pack.inp').write_text(inp)
    with (root/'pack.inp').open() as pack_input:
        result=subprocess.run([str(HERE/'amber-env/bin/packmol')],stdin=pack_input,cwd=root,capture_output=True,text=True,timeout=120)
    (root/'packmol.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0 and (root/'packed.xyz').exists()
    b=read(root/'packed.xyz');assert len(b)==12*args.molecules
    assert np.array_equal(b.numbers,np.tile(a.numbers,args.molecules))
    b.set_cell(np.eye(3)*length);b.pbc=True
    delta=b.positions[:,None,:]-b.positions[None,:,:];delta-=length*np.rint(delta/length)
    distance=np.linalg.norm(delta,axis=-1)
    ids=np.repeat(np.arange(args.molecules),12);mask=ids[:,None]!=ids[None,:]
    nearest=float(distance[mask].min());assert nearest>1.995,nearest
    b.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    write(root/'initial.xyz',b,format='extxyz')
    record=dict(model_sha256=SHA,molecules=args.molecules,atoms=len(b),cell_A=length,
                model_local_cutoff_A=rmax,seed=args.seed,molar_mass_g_mol=molar_mass,
                density_used_only_for_initialization_g_cm3=rho,gas_monomer_fmax_eV_A=fmax,
                periodic_intermolecular_min_distance_A=nearest,
                input_sha256=hashlib.sha256((root/'initial.xyz').read_bytes()).hexdigest(),
                scope='Prepared initial state, no liquid property measured or inferred',
                equilibrium_liquid_qualified=False,cleaning_solution_established=False)
    (root/'prepared.json').write_text(json.dumps(record,indent=2)+'\n');print(record,flush=True)


if __name__=='__main__':main()
