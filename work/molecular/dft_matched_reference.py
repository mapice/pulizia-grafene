"""Matched omegaB97M-V/def2-TZVPD gas reference for learned potentials.

VV10 nonlocal correlation is included by PySCF's functional; do NOT add D3.
Frozen-fragment counterpoise. No solvation entropy or cleaning prediction.
"""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[v] = '1'
os.nice(8)
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
from ase.io import read
from pyscf import dft, gto, lib
lib.num_threads(1)
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--pair',default='HFIP__MethylAcetate')
    p.add_argument('--nlc-grid',type=int,choices=[0,1,2,3],default=0)
    p.add_argument('--initial-grid',type=int,choices=[0,1,2,3],default=None)
    p.add_argument('--monomers-first',action='store_true')
    p.add_argument('--scf-wall-minutes',type=float,default=10.)
    p.add_argument('--resume-progress',action='store_true')
    args=p.parse_args()
    assert 5<=args.scf_wall_minutes<=30
    rec=json.loads((HERE/'results/clusters/binding-refined.json').read_text())['pairs'][args.pair]
    atoms=read(HERE/'results/clusters'/rec['best_geometry'])
    n=next(t['n_acceptor'] for t in rec['trials'] if t['name']+'.refined.xyz'==rec['best_geometry'])
    root=HERE/'results'/f'dft-matched-nlc{args.nlc_grid}';root.mkdir(exist_ok=True)
    energies={}
    sequence=[('AB',range(len(atoms))),('Acp',range(n)),('Bcp',range(n,len(atoms))),
              ('A',range(n)),('B',range(n,len(atoms)))]
    if args.monomers_first:
        sequence=[sequence[3],sequence[4],sequence[0],sequence[1],sequence[2]]
    for tag,real in sequence:
        current=atoms[:n] if tag=='A' else atoms[n:] if tag=='B' else atoms
        real=set(range(len(current))) if tag in ['A','B'] else set(real)
        fingerprint=hashlib.sha256(current.positions.tobytes()+bytes(current.numbers)).hexdigest()
        label=args.pair+'_'+tag;path=root/(label+'.json')
        if path.exists():
            old=json.loads(path.read_text())
            assert old['geometry_sha256']==fingerprint and old['basis']=='def2-tzvpd'
            if old['SCF_converged']:
                energies[tag]=old;continue
            raise RuntimeError('Unconverged point preserved; explicit repair required')
        geometry=[(symbol if i in real else 'ghost-'+symbol,tuple(x))
                  for i,(symbol,x) in enumerate(zip(current.get_chemical_symbols(),current.positions))]
        mol=gto.M(atom=geometry,basis='def2-tzvpd',unit='Angstrom',charge=0,spin=0,
                  verbose=4,max_memory=1800,output=str(root/(label+'.scf.log')))
        mf=dft.RKS(mol).density_fit();mf.xc='wb97m-v'
        mf.grids.level=3;mf.nlcgrids.level=args.nlc_grid;mf.conv_tol=1e-9;mf.max_cycle=100;mf.chkfile=None
        assert mf.do_nlc() and dft.libxc.is_nlc(mf.xc)
        started=time.monotonic()
        def callback(env):
            # An unconverged iterate is restart data, never an accepted energy.
            np.savez(root/(label+'-progress.npz'),density=env['dm'])
            (root/(label+'-progress.json')).write_text(json.dumps(dict(
                SCF_converged=False,geometry_sha256=fingerprint,basis='def2-tzvpd',
                cycle=int(env['cycle']),energy_Ha=float(env['e_tot']),
                elapsed_seconds=time.monotonic()-started),indent=2)+'\n')
            if time.monotonic()-started>args.scf_wall_minutes*60:raise TimeoutError('Bounded SCF wall time')
        mf.callback=callback
        initial=None
        if args.initial_grid is not None:
            cache=HERE/'results'/f'dft-matched-nlc{args.initial_grid}'/(label+'-state.npz')
            previous=json.loads(cache.with_name(label+'.json').read_text())
            assert previous['SCF_converged'] and previous['geometry_sha256']==fingerprint
            initial=np.load(cache)['density']
            assert initial.shape==(mol.nao,mol.nao)
        initial_kind='atomic default' if initial is None else 'previous converged grid'
        if initial is None and args.resume_progress and (root/(label+'-progress.npz')).exists():
            previous=json.loads((root/(label+'-progress.json')).read_text())
            assert previous['geometry_sha256']==fingerprint and previous['basis']=='def2-tzvpd'
            initial=np.load(root/(label+'-progress.npz'))['density'];initial_kind='explicit saved unconverged iterate'
            assert initial.shape==(mol.nao,mol.nao)
        if initial is None and args.monomers_first and tag in ['AB','Acp','Bcp']:
            da=np.load(root/(args.pair+'_A-state.npz'))['density']
            db=np.load(root/(args.pair+'_B-state.npz'))['density']
            assert mol.nao==len(da)+len(db)
            initial=np.zeros((mol.nao,mol.nao))
            if tag in ['AB','Acp']:initial[:len(da),:len(da)]=da
            if tag in ['AB','Bcp']:initial[len(da):,len(da):]=db
            initial_kind='converged isolated-fragment block density'
        e=float(mf.kernel(dm0=initial))
        result=dict(functional='omegaB97M-V with VV10 nonlocal correlation',basis='def2-tzvpd',
                    energy_Ha=e,SCF_converged=bool(mf.converged),no_additional_D3=True,
                    geometry_sha256=fingerprint,basis_functions=mol.nao,real_atoms=len(real),
                    ghost_atoms=len(current)-len(real),grid_level=3,nlc_grid_level=args.nlc_grid,
                    starting_density_grid=args.initial_grid,
                    initial_density_kind=initial_kind,SCF_wall_cap_minutes=args.scf_wall_minutes,
                    elapsed_seconds=time.monotonic()-started)
        path.write_text(json.dumps(result,indent=2)+'\n');print(label,result,flush=True)
        assert mf.converged
        np.savez(root/(label+'-state.npz'),density=mf.make_rdm1(),mo_coeff=mf.mo_coeff,mo_occ=mf.mo_occ)
        energies[tag]=result
    value=lambda k:energies[k]['energy_Ha']
    result=dict(pair=args.pair,geometry=rec['best_geometry'],functional='omegaB97M-V',basis='def2-tzvpd',
                interaction_counterpoise_kJ_mol=(value('AB')-value('Acp')-value('Bcp'))*2625.499638,
                interaction_uncorrected_kJ_mol=(value('AB')-value('A')-value('B'))*2625.499638,
                all_five_SCF_converged=True,energies=energies,nlc_grid_level=args.nlc_grid,
                nonlocal_grid_convergence_not_yet_established=True,
                scope='Matched frozen gas reference, not liquid free energy or polymer removal')
    (root/(args.pair+'_interaction.json')).write_text(json.dumps(result,indent=2)+'\n')
    print('Completed',args.pair,result['interaction_counterpoise_kJ_mol'],flush=True)


if __name__=='__main__':
    main()
