"""Export the same published HFIP model to a second molecular engine.

All harmonic prefactors, units, 1-4 pairs and the unresolved scaling choice
are explicit. Bulk qualification only; no residue-cleaning claim.
"""
from pathlib import Path
import itertools,json,argparse
import numpy as np
from scipy.spatial.transform import Rotation
from ase.io import read
from rdkit import Chem
from published_hfip import MOL,TYPES,NONBONDED,BONDS,ANGLES,canonical,bonded_type,MASSES

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--scale14',type=float,default=.5);parser.add_argument('--side',type=int,default=5);parser.add_argument('--seed',type=int,default=20261006)
    args=parser.parse_args();n=args.side**3
    out=HERE/'results/gromacs'/f'bulk_scale{args.scale14:g}_seed{args.seed}';out.mkdir(parents=True,exist_ok=True)
    top=['; Marchelli2022 SI Tables I-V; kcal/A harmonic prefactors converted explicitly','[ defaults ]',f'1 2 yes {args.scale14:g} {args.scale14:g}','[ atomtypes ]']
    for t,(q,sigma,eps) in NONBONDED.items():
        elem='C' if t in ['CT','CF'] else 'O' if t=='OH' else 'F' if t=='F' else 'H'
        top.append(f'{t:6s} {MASSES[elem]:10.5f} {q:10.6f} A {sigma:.9f} {eps:.9f}')
    top+=['[ moleculetype ]','HFIP 3','[ atoms ]']
    for i,(a,t) in enumerate(zip(MOL.GetAtoms(),TYPES),1):
        top.append(f'{i:3d} {t:5s} 1 HFIP {a.GetSymbol()+str(i):5s} {i:3d} {NONBONDED[t][0]: .7f} {MASSES[a.GetSymbol()]:.5f}')
    top+=['[ bonds ]']
    for b in MOL.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();key=tuple(sorted((bonded_type(TYPES[i]),bonded_type(TYPES[j]))));r,k=BONDS[key]
        top.append(f'{i+1} {j+1} 1 {r:.9f} {k:.8f}')
    top+=['[ angles ]']
    for a in MOL.GetAtoms():
        j=a.GetIdx()
        for i,k in itertools.combinations([x.GetIdx() for x in a.GetNeighbors()],2):
            theta,force=ANGLES[canonical([bonded_type(TYPES[x]) for x in [i,j,k]])]
            top.append(f'{i+1} {j+1} {k+1} 1 {np.rad2deg(theta):.8f} {force:.8f}')
    top+=['[ pairs ]']
    dist=Chem.GetDistanceMatrix(MOL)
    for i,j in itertools.combinations(range(MOL.GetNumAtoms()),2):
        if dist[i,j]==3:top.append(f'{i+1} {j+1} 1')
    top+=['[ system ]','HFIP published-force-field qualification','[ molecules ]',f'HFIP {n}']
    (out/'system.top').write_text('\n'.join(top)+'\n')
    mass=sum(MASSES[a.GetSymbol()] for a in MOL.GetAtoms());box=(n*mass/(602.214076*.8))**(1/3)
    a=read(HERE/'results/clusters/HFIP.refined.xyz');local=a.positions*.1;local-=local.mean(axis=0);rng=np.random.default_rng(args.seed)
    gro=['HFIP model qualification, initial dilute box',f'{n*len(a):5d}'];idx=0
    for m,(i,j,k) in enumerate(itertools.product(range(args.side),repeat=3),1):
        xyz=Rotation.random(random_state=rng).apply(local)+(np.array([i,j,k])+.5)*box/args.side
        for atom,p in zip(MOL.GetAtoms(),xyz):
            idx+=1;name=atom.GetSymbol()+str(atom.GetIdx()+1)
            gro.append(f'{m:5d}{"HFIP":<5s}{name:>5s}{idx:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}')
    gro.append(f'{box:10.5f}{box:10.5f}{box:10.5f}');(out/'initial.gro').write_text('\n'.join(gro)+'\n')
    common='''cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
fourierspacing = 0.12
ewald-rtol = 1e-4
pme-order = 4
constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 6
lincs-iter = 2
'''
    (out/'min.mdp').write_text(common+'''integrator = steep
nsteps = 5000
emtol = 500
emstep = 0.005
nstenergy = 100
''')
    dynamic=common+f'''integrator = md
dt = 0.001
tcoupl = v-rescale
tc-grps = System
tau-t = 0.2
ref-t = 298.15
nstenergy = 1000
nstlog = 1000
nstxout-compressed = 1000
compressed-x-precision = 10000
gen-seed = {args.seed}
'''
    (out/'nvt.mdp').write_text(dynamic+'''nsteps = 50000
pcoupl = no
gen-vel = yes
gen-temp = 298.15
''')
    (out/'npt.mdp').write_text(dynamic+'''nsteps = 1000000
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
gen-vel = no
continuation = yes
''')
    (out/'manifest.json').write_text(json.dumps(dict(source_doi='10.1002/cphc.202100620',source_tables='SI I-V',molecules=n,temperature_K=298.15,initial_density_g_cm3=.8,molar_mass_g_mol=mass,scale14=args.scale14,scale14_verified_from_source=False,scope='Bulk solvent qualification, no PMMA/graphene',compressibility_note='4.5e-5/bar is a chosen barostat tuning parameter, not a measured HFIP material property'),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
