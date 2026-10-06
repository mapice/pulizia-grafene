"""Direct long-distance energy audit for two neutral CF4 molecules.

CF4 is a fluorinated, nonpolar diagnostic with no H bonds. It is not a
cleaning candidate and cannot by itself explain the HFIP equation of state.
Counterpoise wB97M-V/TZVPD uses VV10; no D3 is added to either Hamiltonian.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import hashlib, json, time, subprocess,sys,argparse
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.set_default_dtype(torch.float64)
from ase import Atoms
from ase.io import write,read
from scipy.spatial.transform import Rotation
from mace.calculators import MACECalculator
from pyscf import dft,gto,lib
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products
lib.num_threads(1)

HERE=Path(__file__).resolve().parent
ROOT=HERE/'results/cf4-longrange-audit'
HA_KJ=2625.499638;EV_KJ=96.4853321233


def main():
    global ROOT
    parser=argparse.ArgumentParser();parser.add_argument('--reference-distance',type=float,choices=[8.5,12.],default=8.5)
    args=parser.parse_args()
    if args.reference_distance!=8.5:ROOT=ROOT.with_name(ROOT.name+f'-R{args.reference_distance:.1f}')
    assert not ROOT.exists(),'Preserve completed/interrupted reference'
    ROOT.mkdir(parents=True)
    # TBLite and torch ship different OpenMP runtimes on this Mac. Isolate
    # geometry preparation; never allow duplicate runtimes via KMP flags.
    subprocess.run([sys.executable,str(HERE/'prepare_cf4_geometry.py'),str(ROOT)],check=True,timeout=120)
    geometry_record=json.loads((ROOT/'geometry.json').read_text())
    converged=geometry_record['GFN2_geometry_converged'];assert converged
    a=read(ROOT/'monomer.xyz')
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).double()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')

    def energy(x):
        b=x.copy();b.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
        data=calc._atoms_to_batch(b).to_dict()
        with torch.no_grad():out=model(data,compute_force=False,compute_stress=False,compute_virials=False,training=False)
        return float(out['energy'].item())*calc.energy_units_to_eV

    e0=energy(a);rows=[];reference=None
    for orientation in [0.,90.]:
        b=a.copy();b.positions=b.positions@Rotation.from_euler('z',orientation,degrees=True).as_matrix().T
        em=energy(b)
        for distance in [7.8,8.5,10.,12.,16.,20.]:
            c=b.copy();c.positions[:,0]+=distance;ab=a+c
            nearest=float(np.linalg.norm(a.positions[:,None,:]-c.positions[None,:,:],axis=-1).min())
            assert nearest>float(model.r_max)
            value=(energy(ab)-e0-em)*EV_KJ
            row=dict(carbon_distance_A=distance,orientation_z_degrees=orientation,
                     nearest_cross_atom_A=nearest,all_cross_edges_beyond_local_cutoff=True,
                     interaction_M_kJ_mol=value,rotated_monomer_energy_difference_kJ_mol=(em-e0)*EV_KJ)
            rows.append(row);print(row,flush=True)
            if orientation==0. and distance==args.reference_distance:reference=ab.copy();write(ROOT/'reference-dimer.xyz',reference)
    (ROOT/'learned-curve.json').write_text(json.dumps(dict(model_sha256=SHA,rows=rows,
        local_cutoff_A=float(model.r_max),GFN2_geometry_converged=converged,
        scope='Fixed neutral CF4 gas diagnostic; not an HFIP density or graphene-cleaning prediction'),indent=2)+'\n')
    # Match the training functional on the same frozen geometry. Ghost
    # fragments are calculated independently, not identified by symmetry.
    results=[];old_density={}
    for level in [0,1]:
        energies={}
        for tag,indices in [('A',range(5)),('B',range(5,10)),('AB',range(10)),('Acp',range(5)),('Bcp',range(5,10))]:
            current=(reference[:5] if tag=='A' else reference[5:] if tag=='B' else reference)
            real=set(range(5)) if tag in ['A','B'] else set(indices)
            geometry=[(symbol if i in real else 'ghost-'+symbol,tuple(x))
                      for i,(symbol,x) in enumerate(zip(current.get_chemical_symbols(),current.positions))]
            label=f'{tag}-nlc{level}'
            mol=gto.M(atom=geometry,basis='def2-tzvpd',unit='Angstrom',charge=0,spin=0,
                      verbose=4,max_memory=1600,output=str(ROOT/(label+'.scf.log')))
            mf=dft.RKS(mol).density_fit();mf.xc='wb97m-v';mf.grids.level=3;mf.nlcgrids.level=level
            mf.conv_tol=1e-10;mf.max_cycle=100;mf.chkfile=None
            assert mf.do_nlc() and dft.libxc.is_nlc(mf.xc)
            start=time.monotonic()
            def callback(env):
                np.savez(ROOT/(label+'-progress.npz'),density=env['dm'])
                if time.monotonic()-start>600:raise TimeoutError('SCF own 10-minute cap')
            mf.callback=callback
            initial=old_density.get(tag)
            if initial is None and tag in ['AB','Acp','Bcp']:
                da=old_density['A'];db=old_density['B'];initial=np.zeros((mol.nao,mol.nao))
                if tag in ['AB','Acp']:initial[:len(da),:len(da)]=da
                if tag in ['AB','Bcp']:initial[len(da):,len(da):]=db
            value=float(mf.kernel(dm0=initial));energies[tag]=value
            record=dict(energy_Ha=value,SCF_converged=bool(mf.converged),basis='def2-tzvpd',
                        functional='omegaB97M-V with VV10, no additional D3',nlc_grid_level=level,
                        basis_functions=mol.nao,elapsed_seconds=time.monotonic()-start)
            (ROOT/(label+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(label,record,flush=True)
            assert mf.converged
            old_density[tag]=mf.make_rdm1();np.savez(ROOT/(label+'-state.npz'),density=mf.make_rdm1())
        result=dict(nlc_grid_level=level,energies_Ha=energies,
                    frozen_counterpoise_interaction_kJ_mol=(energies['AB']-energies['Acp']-energies['Bcp'])*HA_KJ,
                    uncorrected_interaction_kJ_mol=(energies['AB']-energies['A']-energies['B'])*HA_KJ,
                    ghost_fragment_symmetry_difference_kJ_mol=(energies['Acp']-energies['Bcp'])*HA_KJ)
        results.append(result);(ROOT/'quantum-progress.json').write_text(json.dumps(results,indent=2)+'\n')
        print(result,flush=True)
    final=dict(model_sha256=SHA,learned_curve=rows,quantum_reference=results,
               quantum_reference_carbon_distance_A=args.reference_distance,quantum_reference_orientation_z_degrees=0.,
               nonlinear_grid_shift_kJ_mol=results[1]['frozen_counterpoise_interaction_kJ_mol']-results[0]['frozen_counterpoise_interaction_kJ_mol'],
               scope='Neutral CF4 frozen gas diagnostic only; no liquid EOS, PMMA/graphene free energy or chemical cleanup inferred',
               no_added_dispersion_to_learned_model=True,cleaning_solution_established=False)
    (ROOT/'complete.json').write_text(json.dumps(final,indent=2)+'\n')


if __name__=='__main__':main()
