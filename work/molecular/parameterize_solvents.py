"""Generate independent GAFF2/AM1-BCC comparator solvent models locally.

AmberTools26 antechamber -> parmchk2 -> tleap -> ParmEd. No guessed
charges or silent missing-parameter defaults; all logs and inputs kept.
Generates a model for validation, not a proof of solvent-cleaning efficacy.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[v]='2'
os.nice(8)
from pathlib import Path
import json,time,subprocess,hashlib
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from tblite.ase import TBLite

HERE=Path(__file__).resolve().parent
AMBER=HERE/'amber-env'
OUT=HERE/'results/gaff2';OUT.mkdir(parents=True,exist_ok=True)
SOLVENTS={'HFIP':('OC(C(F)(F)F)C(F)(F)F','HFI'),
          'DMF':('CN(C)C=O','DMF'),'THF':('C1CCOC1','THF'),
          'AceticAcid':('CC(=O)O','ACE'),'Anisole':('COc1ccccc1','ANI')}
env=os.environ.copy();env['AMBERHOME']=str(AMBER)
env['PATH']=str(AMBER/'bin')+os.pathsep+env['PATH']

def run(stage,cmd,folder):
    start=time.monotonic()
    with (folder/(stage+'.log')).open('w') as f:
        r=subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'guard_run.py'),
                          '--mb','2000','--minutes','10','--']+cmd,
                         cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=650)
    if r.returncode:raise RuntimeError(stage+' failed; inspect '+str(folder/(stage+'.log')))
    return dict(stage=stage,command=cmd,seconds=time.monotonic()-start)

def main():
    records=[]
    for name,(smiles,residue) in SOLVENTS.items():
        folder=OUT/name;folder.mkdir(exist_ok=True)
        if (folder/'manifest.json').exists():
            records.append(json.loads((folder/'manifest.json').read_text()));continue
        mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
        geom=HERE/'results/clusters'/f'{name}.refined.xyz'
        if not geom.exists():
            par=AllChem.ETKDGv3();par.randomSeed=20261006
            if AllChem.EmbedMolecule(mol,par)<0:raise RuntimeError('Embedding failed')
            AllChem.MMFFOptimizeMolecule(mol,maxIters=1000)
            a=Atoms([x.GetSymbol() for x in mol.GetAtoms()],
                    positions=mol.GetConformer().GetPositions())
            a.calc=TBLite(method='GFN2-xTB',verbosity=0,accuracy=1.)
            opt=BFGS(a,logfile=str(folder/'quantum-geometry.log'),maxstep=.10)
            if not opt.run(fmax=.008,steps=300):raise RuntimeError('New molecule geometry unconverged')
            write(geom,a)
        a=read(geom)
        assert mol.GetNumAtoms()==len(a)
        conf=Chem.Conformer(len(a))
        for i,p in enumerate(a.positions):conf.SetAtomPosition(i,tuple(p))
        mol.AddConformer(conf,assignId=True)
        with Chem.SDWriter(str(folder/'input.sdf')) as w:w.write(mol)
        logs=[]
        logs.append(run('antechamber',[str(AMBER/'bin/antechamber'),'-i','input.sdf',
                    '-fi','sdf','-o','molecule.mol2','-fo','mol2','-c','bcc','-nc','0',
                    '-at','gaff2','-rn',residue,'-s','2','-pf','y'],folder))
        # Preserve the generated charges and enforce the known neutral charge
        # with the minimum squared-distance correction. This is not a fit to
        # density/desorption; its size is bounded and reported for each atom.
        raw=folder/'molecule.mol2';lines=raw.read_text().splitlines()
        (folder/'molecule-unadjusted.mol2').write_text(raw.read_text())
        begin=lines.index('@<TRIPOS>ATOM')+1;indices=[];values=[]
        for idx,line in enumerate(lines[begin:],begin):
            if line.startswith('@<TRIPOS>'):break
            if line.strip():indices.append(idx);values.append(float(line.split()[8]))
        excess=sum(values)
        if abs(excess)>.005:raise RuntimeError('Charge error too large for roundoff/neutrality correction')
        correction=-excess/len(values)
        for idx,value in zip(indices,values):
            cols=lines[idx].split();cols[8]=f'{value+correction:.12f}'
            lines[idx]=' '.join(cols)
        raw.write_text('\n'.join(lines)+'\n')
        logs.append(run('parmchk',[str(AMBER/'bin/parmchk2'),'-i','molecule.mol2',
                    '-f','mol2','-s','gaff2','-o','molecule.frcmod'],folder))
        (folder/'leap.in').write_text('source leaprc.gaff2\nloadamberparams molecule.frcmod\n'
                    'M = loadmol2 molecule.mol2\ncheck M\n'
                    'saveamberparm M molecule.prmtop molecule.inpcrd\nquit\n')
        logs.append(run('tleap',[str(AMBER/'bin/tleap'),'-f','leap.in'],folder))
        # ParmEd conversion runs in the isolated Amber environment.
        converter='import parmed as p; s=p.load_file("molecule.prmtop",xyz="molecule.inpcrd"); '+\
                  's.save("molecule.top",overwrite=False); s.save("molecule.gro",overwrite=False); '+\
                  'print("atoms",len(s.atoms),"charge",sum(a.charge for a in s.atoms)); '+\
                  'print("charges",[(a.name,a.type,a.charge) for a in s.atoms])'
        logs.append(run('convert',[str(AMBER/'bin/python'),'-c',converter],folder))
        mol2=(folder/'molecule.mol2').read_text().splitlines()
        j=mol2.index('@<TRIPOS>ATOM');charges=[]
        for line in mol2[j+1:]:
            if line.startswith('@<TRIPOS>'):break
            if line.strip():charges.append(float(line.split()[8]))
        if len(charges)!=len(a) or abs(sum(charges))>5e-5:
            raise RuntimeError('Wrong atom count/net charge')
        frc=(folder/'molecule.frcmod').read_text()
        if 'ATTN' in frc:raise RuntimeError('Unresolved/low-confidence missing parameter flagged ATTN')
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
               for p in folder.iterdir() if p.is_file()}
        record=dict(name=name,smiles=smiles,residue=residue,atoms=len(a),
                    method='GAFF2/AM1-BCC via AmberTools26',charges_e=charges,
                    charge_total_e=sum(charges),generated_unadjusted_charge_total_e=excess,
                    neutrality_correction_per_atom_e=correction,
                    neutrality_method='Minimum L2 charge correction under sum(q)=0; bounded .005e total; original charges kept',
                    logs=logs,file_sha256=files,
                    scope='Comparator force field generated, not qualified for bulk/interfacial desorption.',
                    source='https://ambermd.org/AmberTools.php')
        (folder/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
        records.append(record);print(name,len(a),sum(charges),flush=True)
    (OUT/'solvent-models.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':main()
