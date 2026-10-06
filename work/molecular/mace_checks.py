"""Benchmark an official pinned MACE-OFF23 model against reference clusters.

Only the hash-verified official academic model is deserialized. No silicon
coverage, no graphene-interface qualification, no confidence from branding.
"""
import os
for var in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[var]='2'
os.nice(8)
from pathlib import Path
import json,hashlib,time,resource,sys
import numpy as np
import torch
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.serialization.add_safe_globals([slice])
from mace.calculators import MACECalculator
from ase.io import read

HERE=Path(__file__).resolve().parent
OUT=HERE/'results/mace';OUT.mkdir(parents=True,exist_ok=True)
WEIGHTS=HERE/'sources/MACE-OFF23_small.model'
SHA='165cce4cfec5a34b9c64d4ebf95de15d71106bb584b7291c8470f0749977c46f'
if hashlib.sha256(WEIGHTS.read_bytes()).hexdigest()!=SHA:raise RuntimeError('Unverified model file')
# MACE sets an unsafe default process-wide. Undo it; the one explicit trusted
# model load below retains the unavoidable serialized-module requirement.
os.environ.pop('TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD',None)

def calc(device,dtype):
    torch.set_default_dtype(torch.float32 if dtype=='float32' else torch.float64)
    model=torch.load(WEIGHTS,map_location='cpu',weights_only=False)
    # Constructor moves model to device before converting dtype. Metal has
    # no float64, so convert on CPU first rather than changing global policy.
    if dtype=='float32':model=model.float()
    return MACECalculator(models=model,device=device,default_dtype=dtype)

def evaluation(a,c):
    a=a.copy();a.calc=c
    start=time.monotonic();e=float(a.get_potential_energy());f=a.get_forces()
    if c.device.type=='mps':torch.mps.synchronize()
    return e,f,time.monotonic()-start

def main():
    data=json.loads((HERE/'results/clusters/binding-refined.json').read_text())
    c=calc('cpu','float64');start=time.monotonic();rows=[]
    for key,r in data['pairs'].items():
        a=read(HERE/'results/clusters'/r['best_geometry'])
        n=next(t['n_acceptor'] for t in r['trials'] if t['name']+'.refined.xyz'==r['best_geometry'])
        e,force,seconds=evaluation(a,c);ea,_,_=evaluation(a[:n],c);eb,_,_=evaluation(a[n:],c)
        row=dict(pair=key,interaction_frozen_kJ_mol=(e-ea-eb)*96.4853321233,max_force_at_GFN_geometry_eV_A=float(np.linalg.norm(force,axis=1).max()),elapsed_complex_seconds=seconds)
        rows.append(row);print(json.dumps(row),flush=True)
    a=read(HERE/'results/clusters'/data['pairs']['HFIP__MethylAcetate']['best_geometry'])
    modes=[('cpu','float64'),('cpu','float32')]
    if torch.backends.mps.is_available():modes.append(('mps','float32'))
    benchmarks=[];reference=None
    for device,dtype in modes:
        try:
            cc=c if device=='cpu' and dtype=='float64' else calc(device,dtype)
            e,f,t=evaluation(a,cc)
            if reference is None:reference=(e,f)
            times=[]
            for k in range(3):
                probe=a.copy();probe.positions[0,0]+=1e-4*(k+1)
                _,_,t=evaluation(probe,cc);times.append(t)
            row=dict(device=device,dtype=dtype,median_seconds=float(np.median(times)),energy_difference_from_CPU64_eV=e-reference[0],forces_RMS_difference_eV_A=float(np.sqrt(np.mean((f-reference[1])**2))))
        except Exception as exc:
            import traceback
            row=dict(device=device,dtype=dtype,error=str(exc),traceback=traceback.format_exc())
        benchmarks.append(row);print(json.dumps(row),flush=True)
    record=dict(model='MACE-OFF23-small',model_sha256=SHA,torch_version=torch.__version__,model_elements=[int(x) for x in c.models[0].atomic_numbers],scope='Gas-phase proxy checks. Model lacks Si; extended graphene not qualified. Training reference differs from PBE0-D3(BJ); differences are diagnostics, not direct universal errors.',interactions=rows,benchmarks=benchmarks,elapsed_seconds=time.monotonic()-start,peak_RSS_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024))
    (OUT/'cluster-checks.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Finished',record['elapsed_seconds'],record['peak_RSS_MB'],flush=True)

if __name__=='__main__':main()
