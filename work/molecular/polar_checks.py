"""Actual MACE-POLAR-1-S checks against molecular proxies and electric field.

Official release digest checked before loading. Molecular training does not
qualify extended graphene or silica. Differences from PBE0 are diagnostics
because the training functional is omegaB97M-V.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='2'
os.nice(8)
from pathlib import Path
import json,time,hashlib,resource,sys
import numpy as np
import torch
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.serialization.add_safe_globals([slice])
from mace.calculators import MACECalculator
from ase.io import read

HERE=Path(__file__).resolve().parent;OUT=HERE/'results/polar';OUT.mkdir(parents=True,exist_ok=True)
FILE=HERE/'sources/MACE-POLAR-1-S.model'
SHA='e4495612037b3b3312633182882a38a694ecac9ea0be2b9889ac0b2a84a99510'
assert hashlib.sha256(FILE.read_bytes()).hexdigest()==SHA
os.environ.pop('TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD',None)
model=torch.load(FILE,map_location='cpu',weights_only=False)
calc=MACECalculator(models=model,device='cpu',default_dtype='float64',model_type='PolarMACE',pbc_handling='realspace')
# graph_longrange0.4.0 matches the released checkpoint/caller API. v0.4.3/4.4
# changed data layout and are not interchangeable; no physics adapter is used.

def energy(a,field=None):
    a=a.copy();a.info.update(charge=0,spin=1,external_field=[0.,0.,0.] if field is None else field);a.calc=calc
    start=time.monotonic();e=float(a.get_potential_energy());forces=a.get_forces()
    r=dict(energy_eV=e,elapsed_seconds=time.monotonic()-start,max_force_eV_A=float(np.linalg.norm(forces,axis=1).max()))
    for k in ['dipole','charges','density_coefficients']:
        if k in calc.results:r[k]=np.array(calc.results[k]).tolist()
    if not np.isfinite(e) or not np.isfinite(forces).all():raise RuntimeError('Nonfinite ML result')
    return r

def main():
    d=json.loads((HERE/'results/clusters/binding-refined.json').read_text());rows=[];start=time.monotonic()
    for key,r in d['pairs'].items():
        a=read(HERE/'results/clusters'/r['best_geometry'])
        n=next(t['n_acceptor'] for t in r['trials'] if t['name']+'.refined.xyz'==r['best_geometry'])
        ab=energy(a);aa=energy(a[:n]);bb=energy(a[n:])
        row=dict(pair=key,interaction_frozen_kJ_mol=(ab['energy_eV']-aa['energy_eV']-bb['energy_eV'])*96.4853321233,complex=ab,acceptor=aa,donor=bb)
        rows.append(row);print(key,row['interaction_frozen_kJ_mol'],ab['elapsed_seconds'],flush=True)
        (OUT/'clusters-partial.json').write_text(json.dumps(rows,indent=2)+'\n')
    a=read(HERE/'results/clusters/HFIP.refined.xyz');fields=[]
    for f in [0.,-.002,.002]:fields.append(dict(field_V_A=[0.,0.,f],**energy(a,[0.,0.,f])))
    # Field response is gas-phase at fixed nuclei, not bulk permittivity.
    response=-(fields[1]['energy_eV']+fields[2]['energy_eV']-2*fields[0]['energy_eV'])/(.002**2)
    record=dict(model='MACE-POLAR-1-S',model_sha256=SHA,elements=[int(x) for x in calc.models[0].atomic_numbers],interactions=rows,HFIP_field_checks=fields,HFIP_fixed_nuclei_energy_curvature_eA2_V=response,scope='Gas electronic proxies and fixed-nuclei field response; neither liquid dielectric nor residue cleaning. Molecular training does not qualify extended graphene/SiO2.',elapsed_seconds=time.monotonic()-start,peak_RSS_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024))
    (OUT/'checks.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Finished',record['elapsed_seconds'],record['peak_RSS_MB'],flush=True)

if __name__=='__main__':main()
