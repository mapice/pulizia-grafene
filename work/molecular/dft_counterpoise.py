"""DFT electronic interaction checks at frozen, refined GFN2 geometries.

PBE0-D3(BJ), density fitting. Boys-Bernardi counterpoise corrects basis
superposition error. Dispersion uses only REAL atoms in each fragment.
This is neither liquid free energy nor a prediction of residue removal.
"""
import os
for var in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[var]='2'
os.nice(8)
from pathlib import Path
import json,time,resource,sys,argparse,hashlib
import numpy as np
from ase.io import read
from pyscf import gto,dft,lib
from dftd3.interface import DispersionModel,RationalDampingParam

lib.num_threads(2)
HERE=Path(__file__).resolve().parent
CLUST=HERE/'results/clusters'
OUT=HERE/'results/dft'
OUT.mkdir(parents=True,exist_ok=True)
HARTREE_KJMOL=2625.499638
BOHR_A=.529177210903

def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)

def energy(atoms,real_indices,label,args):
    path=OUT/f'{label}.json'
    xyz=atoms.positions
    geometry_hash=hashlib.sha256(xyz.tobytes()+bytes(atoms.numbers)).hexdigest()
    if path.exists():
        old=json.loads(path.read_text())
        if old['geometry_hash']==geometry_hash and old['basis']==args.basis and old['xc']==args.xc and old['converged']:
            return old
    real=set(real_indices)
    atom_input=[(symbol if i in real else 'ghost-'+symbol,tuple(p)) for i,(symbol,p) in enumerate(zip(atoms.get_chemical_symbols(),xyz))]
    mol=gto.M(atom=atom_input,basis=args.basis,unit='Angstrom',charge=0,spin=0,verbose=4,max_memory=1800,output=str(OUT/f'{label}.scf.log'))
    mf=dft.RKS(mol).density_fit()
    mf.xc=args.xc
    mf.grids.level=3
    mf.conv_tol=1e-9
    mf.max_cycle=120
    mf.chkfile=None
    start=time.monotonic()
    def callback(env):
        if rss()>3000:raise MemoryError('3 GB RSS limit exceeded')
        if time.monotonic()-start>1800:raise TimeoutError('30 min single-point cap')
    mf.callback=callback
    e=float(mf.kernel())
    real_atoms=atoms[list(real_indices)]
    model=DispersionModel(real_atoms.numbers,real_atoms.positions/BOHR_A)
    correction=float(model.get_dispersion(RationalDampingParam(method=args.xc),grad=False)['energy'])
    r=dict(label=label,basis=args.basis,xc=args.xc,converged=bool(mf.converged),electronic_Ha=e,dispersion_Ha=correction,total_Ha=e+correction,basis_functions=mol.nao,real_atoms=len(real),ghost_atoms=len(atoms)-len(real),geometry_hash=geometry_hash,elapsed_seconds=time.monotonic()-start,peak_rss_MB=rss(),grid_level=3,conv_tol_Ha=1e-9)
    path.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r),flush=True)
    if not mf.converged:raise RuntimeError('SCF did not converge '+label)
    return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--pair',default='HFIP__MethylAcetate')
    p.add_argument('--basis',default='def2-svp')
    p.add_argument('--xc',default='pbe0')
    args=p.parse_args()
    r=json.loads((CLUST/'binding-refined.json').read_text())['pairs'][args.pair]
    atoms=read(CLUST/r['best_geometry'])
    n=next(t['n_acceptor'] for t in r['trials'] if t['name']+'.refined.xyz'==r['best_geometry'])
    tag=f'{args.pair}_{args.xc}_{args.basis}'
    ab=energy(atoms,range(len(atoms)),tag+'_AB',args)
    a_cp=energy(atoms,range(n),tag+'_Acp',args)
    b_cp=energy(atoms,range(n,len(atoms)),tag+'_Bcp',args)
    a=energy(atoms[:n],range(n),tag+'_A',args)
    b=energy(atoms[n:],range(len(atoms)-n),tag+'_B',args)
    ecp=ab['total_Ha']-a_cp['total_Ha']-b_cp['total_Ha']
    eraw=ab['total_Ha']-a['total_Ha']-b['total_Ha']
    record=dict(pair=args.pair,method=args.xc+'-D3(BJ)',basis=args.basis,geometry=r['best_geometry'],geometry_method='GFN2-xTB refined; not DFT optimized',interaction_counterpoise_kJ_mol=ecp*HARTREE_KJMOL,interaction_raw_kJ_mol=eraw*HARTREE_KJMOL,BSSE_correction_kJ_mol=(ecp-eraw)*HARTREE_KJMOL,scope='Vacuum frozen-fragment electronic interaction energies. No liquid entropy/solvation, no cleaning efficacy.',all_converged=all(x['converged'] for x in [ab,a_cp,b_cp,a,b]),peak_rss_MB=rss(),energies=dict(AB=ab,Acp=a_cp,Bcp=b_cp,A=a,B=b))
    (OUT/f'{tag}_interaction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)

if __name__=='__main__':main()
