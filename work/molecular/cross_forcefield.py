"""Cross-check published intermolecular LJ/Coulomb interactions.

Exact published HFIP/acetone atom parameters; an explicitly labelled capped
PMMA local proxy. No force-field fitting to make the chosen solvent win.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='2'
os.nice(8)
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from ase.io import read
from published_hfip import NONBONDED,TYPES

HERE=Path(__file__).resolve().parent
COULOMB=138.93545764438198 # kJ nm mol^-1 e^-2

def params(name):
    if name=='HFIP':return np.array([NONBONDED[t] for t in TYPES]),'published SI Tables I and V'
    if name=='Acetone':
        mol=Chem.AddHs(Chem.MolFromSmiles('CC(=O)C'))
        values={'C':(.47,.375,.105*4.184),'CT':(-.18,.35,.066*4.184),'HC':(.06,.25,.03*4.184),'O':(-.47,.296,.210*4.184)}
        tags=['HC' if a.GetSymbol()=='H' else 'O' if a.GetSymbol()=='O' else 'C' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in a.GetBonds()) else 'CT' for a in mol.GetAtoms()]
        return np.array([values[t] for t in tags]),'published SI Tables I and V'
    if name=='MethylPivalate':
        mol=Chem.AddHs(Chem.MolFromSmiles('CC(C)(C)C(=O)OC'))
        row=[];caps=[]
        for a in mol.GetAtoms():
            i=a.GetIdx();z=a.GetSymbol()
            if z=='O':r=(-.5939,.2959922,.878640) if any(b.GetBondType()==Chem.BondType.DOUBLE for b in a.GetBonds()) else (-.4617,.3000012,.711280)
            elif z=='H':r=(0.,.2471353,.065689) if any(n.GetSymbol()=='O' for n in a.GetNeighbors()[0].GetNeighbors()) else (0.,.2649533,.065689)
            elif any(b.GetBondType()==Chem.BondType.DOUBLE for b in a.GetBonds()):r=(.7464,.339967,.359824)
            elif any(n.GetSymbol()=='O' for n in a.GetNeighbors()):r=(.2801,.339967,.457730)
            elif sum(n.GetSymbol()=='C' for n in a.GetNeighbors())==4:r=(.0189,.339967,.457730)
            else:r=(.0051,.339967,.457730);caps.append(i)
            row.append(r)
        row=np.array(row);excess=float(row[:,0].sum())
        row[caps,0]-=excess/len(caps)
        return row,'PMMA2021 Table S2 types applied to capped proxy; net +0.0051e distributed over three caps, an unvalidated boundary correction'
    raise ValueError(name)

def interaction(xa,pa,xb,pb):
    r=np.linalg.norm(xa[:,None,:]-xb[None,:,:],axis=2)*.1
    sigma=(pa[:,None,1]+pb[None,:,1])/2
    eps=np.sqrt(pa[:,None,2]*pb[None,:,2])
    sr6=(sigma/r)**6
    lj=float(np.sum(4*eps*(sr6**2-sr6)))
    el=float(np.sum(COULOMB*pa[:,None,0]*pb[None,:,0]/r))
    return dict(LJ_kJ_mol=lj,Coulomb_kJ_mol=el,total_kJ_mol=lj+el)

def main():
    d=json.loads((HERE/'results/clusters/binding-refined.json').read_text());rows=[]
    for key in ['HFIP__HFIP','HFIP__Acetone','HFIP__MethylPivalate']:
        r=d['pairs'][key];a=read(HERE/'results/clusters'/r['best_geometry'])
        t=next(t for t in r['trials'] if t['name']+'.refined.xyz'==r['best_geometry']);n=t['n_acceptor']
        pa,ma=params(r['acceptor']);pb,mb=params(r['donor'])
        assert len(pa)==n and len(pb)==len(a)-n
        assert abs(pa[:,0].sum())<1e-10 and abs(pb[:,0].sum())<1e-10
        base=interaction(a.positions[:n],pa,a.positions[n:],pb)
        direction=a.positions[t['donor_H']]-a.positions[t['acceptor_O']];direction/=np.linalg.norm(direction)
        curve=[]
        for delta in np.linspace(-.25,2.0,46):
            curve.append(dict(translation_A=float(delta),**interaction(a.positions[:n],pa,a.positions[n:]+delta*direction,pb)))
        rows.append(dict(pair=key,geometry=r['best_geometry'],baseline=base,acceptor_parameter_scope=ma,donor_parameter_scope=mb,separation_curve=curve))
        print(key,base,flush=True)
    (HERE/'results/cross-forcefield.json').write_text(json.dumps(dict(method='Published Lorentz-Berthelot LJ plus Coulomb at frozen GFN2 geometries',scope='Intermolecular model audit, not solvation/desorption/cleaning prediction',rows=rows),indent=2)+'\n')

if __name__=='__main__':main()
