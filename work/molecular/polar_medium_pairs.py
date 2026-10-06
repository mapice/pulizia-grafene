"""MACE-POLAR-1-M: identical frozen gas complexes used by the other models."""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import json,hashlib,time
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
from ase.io import read
from mace.calculators import MACECalculator
HERE=Path(__file__).resolve().parent
SHA='fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'
path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False)
calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')

def energy(a):
    b=a.copy();b.positions-=b.positions.mean(0)
    b.info.update(charge=0,spin=1,external_field=[0.,0.,0.]);b.calc=calc
    value=float(b.get_potential_energy());assert np.isfinite(value)
    return value

started=time.monotonic();rows=[]
data=json.loads((HERE/'results/clusters/binding-refined.json').read_text())['pairs']
for key,r in data.items():
    a=read(HERE/'results/clusters'/r['best_geometry'])
    n=next(t['n_acceptor'] for t in r['trials'] if t['name']+'.refined.xyz'==r['best_geometry'])
    e=(energy(a)-energy(a[:n])-energy(a[n:]))*96.4853321233
    row=dict(pair=key,geometry=r['best_geometry'],interaction_frozen_kJ_mol=e)
    rows.append(row);print(row,flush=True)
root=HERE/'results/polar-size-audit'
(root/'M-pairs.json').write_text(json.dumps(dict(
    model='MACE-POLAR-1-M',sha256=SHA,pairs=rows,elapsed_seconds=time.monotonic()-started,
    scope='Gas electronic proxies, not a liquid free-energy or cleaning ranking'),indent=2)+'\n')
