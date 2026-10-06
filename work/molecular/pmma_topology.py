"""PMMA oligomer model transcribed from Behbahani/Harmandaris 2021 SI S2.

Every bond/angle/dihedral must match an explicit table entry: missing terms
raise an exception. Neutral end caps are a documented local model choice,
not published parameters or a substitute for a 950K-chain calculation.
"""
from pathlib import Path
import itertools,json,argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import write

HERE=Path(__file__).resolve().parent
def canon(t):return min(tuple(t),tuple(t)[::-1])
NB={'CX':(.0051,.339967,.457730,12.011),'CA':(.0189,.339967,.457730,12.011),'CE':(.2801,.339967,.457730,12.011),'C':(.7464,.339967,.359824,12.011),'O':(-.5939,.2959922,.878640,15.9994),'OS':(-.4617,.3000012,.711280,15.9994),'HC':(0.,.2649533,.065689,1.008),'H1':(0.,.2471353,.065689,1.008)}
BOND={canon(k):(r,f) for k,r,f in [(('CX','CX'),.152,259408.),(('HC','CX'),.1095,301248.),(('CX','C'),.150,276144.),(('C','O'),.1204,835963.2),(('C','OS'),.1343,344175.9),(('OS','CX'),.141,267776.),(('CX','H1'),.1092,343088.)]}
ANGLE={canon(k):(theta,force) for k,theta,force in [(('HC','CX','HC'),109.5,334.7),(('HC','CX','CX'),112.6,376.6),(('CX','CX','CX'),113.5,376.6),(('CX','CX','C'),111.5,334.7),(('CX','C','O'),125.4,493.7),(('CX','C','OS'),111.,418.4),(('C','OS','CX'),114.,423.4),(('O','C','OS'),122.5,861.9),(('OS','CX','H1'),110.,502.1),(('H1','CX','H1'),109.5,376.6)]}
TORSION={canon(k):coef for k,coef in [(('CX','CX','CX','CX'),(1.8828,0,0)),(('CX','CX','CX','HC'),(0,0,.4184)),(('CX','CX','CX','C'),(0,0,0)),(('HC','CX','CX','C'),(0,0,-.12552)),(('CX','CX','C','O'),(0,-1.7154,-1.2134)),(('CX','CX','C','OS'),(0,-1.8828,-1.3807)),(('O','C','OS','CX'),(0,-17.61464,-3.05432)),(('H1','CX','OS','C'),(0,0,.6276)),(('CX','C','OS','CX'),(-2.7196,-11.506,0))]}
def category(t):return 'CX' if t in ['CA','CE'] else t

def molecule(n,seed=20261006):
    rng=np.random.default_rng(seed)
    symbols=['[C@]' if b else '[C@@]' for b in rng.integers(0,2,n)]
    smiles='C'+''.join(c+'(C)(C(=O)OC)C' for c in symbols)
    heavy=Chem.MolFromSmiles(smiles);last=heavy.GetNumAtoms()-1
    mol=Chem.AddHs(heavy)
    tags=[]
    for a in mol.GetAtoms():
        z=a.GetSymbol()
        if z=='O':t='O' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in a.GetBonds()) else 'OS'
        elif z=='H':t='H1' if any(x.GetSymbol()=='O' for x in a.GetNeighbors()[0].GetNeighbors()) else 'HC'
        elif any(b.GetBondType()==Chem.BondType.DOUBLE for b in a.GetBonds()):t='C'
        elif any(x.GetSymbol()=='O' for x in a.GetNeighbors()):t='CE'
        elif a.GetDegree()==4 and all(x.GetSymbol()=='C' for x in a.GetNeighbors()):t='CA'
        else:t='CX'
        tags.append(t)
    charges=np.array([NB[t][0] for t in tags]);excess=float(charges.sum())
    assert np.isclose(excess,.0051,atol=1e-12)
    charges[[0,last]]-=excess/2
    assert abs(charges.sum())<1e-12
    params=AllChem.ETKDGv3();params.randomSeed=seed
    ids=AllChem.EmbedMultipleConfs(mol,numConfs=4,params=params)
    if not ids:raise RuntimeError('Oligomer embedding failed')
    opt=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000)
    cid=int(ids[min(range(len(ids)),key=lambda i:opt[i][1])])
    xyz=np.array(mol.GetConformer(cid).GetPositions())
    return mol,tags,charges,xyz,smiles,[0,last]

def export(n,seed=20261006):
    out=HERE/'results/pmma'/f'n{n}_seed{seed}';out.mkdir(parents=True,exist_ok=True)
    mol,tags,charges,xyz,smiles,caps=molecule(n,seed)
    top=['; Source10.3390/polym13050830 SI TableS2; excluded 1-2/1-3, unscaled 1-4 assumption','[ moleculetype ]','PMMA 2','[ atoms ]']
    atomtypes=['[ atomtypes ]']
    for t,(q,sigma,eps,mass) in NB.items():atomtypes.append(f'P{t} {mass:.6f} {q:.7f} A {sigma:.9f} {eps:.9f}')
    for i,(t,q) in enumerate(zip(tags,charges),1):top.append(f'{i} P{t} 1 PMMA {t+str(i):6s} {i} {q:.9f} {NB[t][3]:.6f}')
    top+=['[ bonds ]'];bond_records=[]
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();key=canon([category(tags[x]) for x in [i,j]]);r,k=BOND[key]
        top.append(f'{i+1} {j+1} 1 {r:.9f} {k:.9f}');bond_records.append((i,j,r,k))
    top+=['[ angles ]'];angle_records=[]
    for a in mol.GetAtoms():
        j=a.GetIdx()
        for i,k in itertools.combinations([x.GetIdx() for x in a.GetNeighbors()],2):
            key=canon([category(tags[x]) for x in [i,j,k]]);theta,f=ANGLE[key]
            top.append(f'{i+1} {j+1} {k+1} 1 {theta:.9f} {f:.9f}');angle_records.append((i,j,k,theta,f))
    top+=['[ dihedrals ]'];torsion_records=[]
    for b in mol.GetBonds():
        j,k=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        for i in [x.GetIdx() for x in mol.GetAtomWithIdx(j).GetNeighbors() if x.GetIdx()!=k]:
            for l in [x.GetIdx() for x in mol.GetAtomWithIdx(k).GetNeighbors() if x.GetIdx()!=j]:
                key=canon([category(tags[x]) for x in [i,j,k,l]])
                coef=TORSION[key] # Unknown entries are never set to zero silently.
                for mult,f in enumerate(coef,1):
                    if f:
                        top.append(f'{i+1} {j+1} {k+1} {l+1} 9 0.0 {f:.9f} {mult}')
                        torsion_records.append((i,j,k,l,mult,f))
    (out/'pmma.itp').write_text('\n'.join(top)+'\n');(out/'atomtypes.itp').write_text('\n'.join(atomtypes)+'\n')
    write(out/'initial.xyz',Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=xyz))
    manifest=dict(repeat_units=n,atoms=mol.GetNumAtoms(),smiles=smiles,seed=seed,tacticity='seeded mixed stereocenters, representative atactic oligomer only',charge_total_e=float(charges.sum()),end_cap_atoms_zero_based=caps,end_correction_e_each=-.00255,end_correction_validated=False,source_doi='10.3390/polym13050830',source_SI_table='S2',bond_count=len(bond_records),angle_count=len(angle_records),nonzero_torsion_count=len(torsion_records),nonbonded_exclusion='1-2,1-3 excluded; 1-4 fully unscaled assumption pending sensitivity audit',scope='Oligomer model built. Not equilibrated; not PMMA950K, no cleaning prediction.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(manifest)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--seed',type=int,default=20261006);a=p.parse_args();export(a.n,a.seed)
