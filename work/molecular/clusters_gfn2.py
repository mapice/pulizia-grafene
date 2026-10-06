"""Quantum tight-binding molecular checks, not solvent cleaning predictions.

GFN2-xTB electronic binding/interaction energies in VACUUM. Polymer proxies
are small molecules, not PMMA950K or PPC chains. Every quantity is labelled.
"""
import os
for var in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
 os.environ[var]='2'
os.nice(8)
import numpy as np
from scipy.spatial.transform import Rotation
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write,read
from tblite.ase import TBLite
from pathlib import Path
import json,time,resource,sys,argparse

HERE=Path(__file__).resolve().parent
OUT=HERE/'results/clusters'
OUT.mkdir(parents=True,exist_ok=True)
EV_KJMOL=96.4853321233
MOLECULES={
 'HFIP':'OC(C(F)(F)F)C(F)(F)F',
 'IPA':'CC(O)C',
 'AceticAcid':'CC(=O)O',
 'TFE':'OCC(F)(F)F',
 'MeOH':'CO',
 'Water':'O',
 'Acetone':'CC(=O)C',
 'DMF':'CN(C)C=O',
 'THF':'C1CCOC1',
 'MethylAcetate':'CC(=O)OC',
 'MethylPivalate':'CC(C)(C)C(=O)OC',
 'DimethylCarbonate':'COC(=O)OC',
}
PAIRS=[('HFIP',a) for a in ['HFIP','MeOH','Acetone','DMF','THF','MethylAcetate','MethylPivalate','DimethylCarbonate']]
PAIRS += [(d,'MethylAcetate') for d in ['IPA','AceticAcid','TFE','MeOH','Water']]

def calculator():
 return TBLite(method='GFN2-xTB',verbosity=0,accuracy=1.0,max_iterations=250)

def rss_mb():
 return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)

def optimize(atoms,name):
 atoms.calc=calculator()
 opt=BFGS(atoms,logfile=str(OUT/f'{name}.opt.log'),trajectory=str(OUT/f'{name}.traj'),maxstep=.15)
 converged=bool(opt.run(fmax=.025,steps=100))
 energy=float(atoms.get_potential_energy())
 force=float(np.linalg.norm(atoms.get_forces(),axis=1).max())
 write(OUT/f'{name}.xyz',atoms)
 if not np.isfinite(energy) or not np.isfinite(atoms.positions).all(): raise RuntimeError('Non-finite geometry/energy')
 if rss_mb()>3000: raise MemoryError('3 GB RSS limit exceeded')
 return energy,converged,force

def molecule(name):
 mol=Chem.AddHs(Chem.MolFromSmiles(MOLECULES[name]))
 parameters=AllChem.ETKDGv3()
 parameters.randomSeed=20261006
 ids=AllChem.EmbedMultipleConfs(mol,numConfs=6,params=parameters)
 if not ids: raise RuntimeError('Embedding failed')
 energies=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=500,mmffVariant='MMFF94s')
 cid=int(ids[min(range(len(ids)),key=lambda i:energies[i][1])])
 positions=np.array(mol.GetConformer(cid).GetPositions())
 atoms=Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=positions)
 return mol,atoms

def donor_indices(mol):
 for atom in mol.GetAtoms():
  if atom.GetSymbol()=='O':
   hs=[n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol()=='H']
   if hs:return atom.GetIdx(),hs[0]
 raise ValueError('No OH donor')

def acceptor_indices(mol):
 for bond in mol.GetBonds():
  atoms=[bond.GetBeginAtom(),bond.GetEndAtom()]
  if bond.GetBondType()==Chem.BondType.DOUBLE and {a.GetSymbol() for a in atoms}=={'C','O'}:
   return next(a.GetIdx() for a in atoms if a.GetSymbol()=='O')
 return next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=='O')

def starting_complex(dmol,donor,amol,acceptor,angle):
 ao=acceptor_indices(amol)
 a=acceptor.copy()
 neighbors=[n.GetIdx() for n in amol.GetAtomWithIdx(ao).GetNeighbors() if n.GetSymbol()!='H']
 if not neighbors:neighbors=[n.GetIdx() for n in amol.GetAtomWithIdx(ao).GetNeighbors()]
 direction=a.positions[ao]-a.positions[neighbors].mean(axis=0)
 direction/=np.linalg.norm(direction)
 a.positions-=a.positions[ao].copy()
 # Align acceptor's outward direction with +x. Row-vector convention.
 rot,_=Rotation.align_vectors(np.array([[1.,0,0]]),direction.reshape(1,3))
 a.positions=rot.apply(a.positions)
 do,dh=donor_indices(dmol)
 d=donor.copy()
 axis=d.positions[dh]-d.positions[do]
 axis/=np.linalg.norm(axis)
 rot,_=Rotation.align_vectors(np.array([[-1.,0,0]]),axis.reshape(1,3))
 d.positions-=d.positions[dh].copy()
 d.positions=rot.apply(d.positions)
 d.positions=Rotation.from_rotvec(np.array([angle,0,0])).apply(d.positions)
 d.positions+=np.array([1.75,0,0])
 return a+d,dict(acceptor_O=ao,donor_O=len(a)+do,donor_H=len(a)+dh,n_acceptor=len(a))

def singlepoint(atoms):
 atoms.calc=calculator()
 return float(atoms.get_potential_energy())

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--pair',default=None,help='donor:acceptor, otherwise all')
 args=parser.parse_args()
 pairs=[tuple(args.pair.split(':'))] if args.pair else PAIRS
 monomers_path=OUT/'monomers.json'
 monomers=json.loads(monomers_path.read_text()) if monomers_path.exists() else {}
 objects={}
 start=time.monotonic()
 for name in dict.fromkeys(x for pair in pairs for x in pair):
  mol,atoms=molecule(name)
  if name in monomers and (OUT/f'{name}.xyz').exists():
   atoms=read(OUT/f'{name}.xyz')
  else:
   e,c,f=optimize(atoms,name)
   monomers[name]=dict(smiles=MOLECULES[name],atoms=len(atoms),energy_eV=e,converged=c,max_force_eV_A=f)
   monomers_path.write_text(json.dumps(monomers,indent=2)+'\n')
  objects[name]=(mol,atoms)
  print(json.dumps(dict(monomer=name,**monomers[name])),flush=True)
 result_path=OUT/'binding.json'
 data=json.loads(result_path.read_text()) if result_path.exists() else {'method':'GFN2-xTB','scope':'Vacuum electronic energies of molecular proxies. Not liquid free energies; not cleaning/solvent ranking.','pairs':{}}
 for dname,aname in pairs:
  key=dname+'__'+aname
  if key in data['pairs'] and data['pairs'][key].get('finished'):continue
  dmol,donor=objects[dname]; amol,acceptor=objects[aname]
  trials=[]
  for index,angle in enumerate([0.,2*np.pi/3,4*np.pi/3]):
   name=key+f'_{index}'
   try:
    complex_,indices=starting_complex(dmol,donor,amol,acceptor,angle)
    e,c,f=optimize(complex_,name)
    trials.append(dict(name=name,energy_eV=e,converged=c,max_force_eV_A=f,**indices))
   except Exception as exc:
    trials.append(dict(name=name,error=str(exc)))
   print(json.dumps(dict(pair=key,trial=trials[-1],elapsed_seconds=time.monotonic()-start,peak_rss_mb=rss_mb())),flush=True)
  valid=[x for x in trials if x.get('converged')]
  if not valid:
   data['pairs'][key]=dict(finished=False,trials=trials,error='No converged start')
   result_path.write_text(json.dumps(data,indent=2)+'\n');continue
  best=min(valid,key=lambda x:x['energy_eV'])
  ab=read(OUT/f"{best['name']}.xyz")
  n=best['n_acceptor']
  ea=singlepoint(ab[:n]); ed=singlepoint(ab[n:])
  eint=best['energy_eV']-ea-ed
  eb=best['energy_eV']-monomers[aname]['energy_eV']-monomers[dname]['energy_eV']
  oh=float(np.linalg.norm(ab.positions[best['donor_O']]-ab.positions[best['donor_H']]))
  h_o=float(np.linalg.norm(ab.positions[best['acceptor_O']]-ab.positions[best['donor_H']]))
  record=dict(finished=True,donor=dname,acceptor=aname,best_geometry=best['name']+'.xyz',binding_electronic_kJ_mol=eb*EV_KJMOL,interaction_frozen_fragments_kJ_mol=eint*EV_KJMOL,deformation_kJ_mol=(eb-eint)*EV_KJMOL,OH_A=oh,H_acceptor_A=h_o,trials=trials,peak_rss_mb=rss_mb(),elapsed_seconds=time.monotonic()-start)
  data['pairs'][key]=record
  result_path.write_text(json.dumps(data,indent=2)+'\n')
  print(json.dumps(dict(result=key,**record)),flush=True)
 print('Finished elapsed_seconds',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
