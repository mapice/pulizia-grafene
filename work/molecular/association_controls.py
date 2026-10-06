"""Known canonical distribution for an angular association volume.

Flat and piecewise-constant reference potentials have exact populations.
Omission of the Hastings factor is an explicitly labelled negative control.
This is NOT a calculation of the HFIP fluid or a cleaning prediction.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import json,time
import numpy as np
from association_proposal import propose,eligible_acceptors

HERE=Path(__file__).resolve().parent


def run(binding,corrected,seed,steps=100000):
    rng=np.random.default_rng(seed);L=15.;rlo=2.6;rhi=3.3;cmin=np.cos(np.pi/6);bias=.5
    # Labelled rigid trihedra test translational/angular measure only.
    shape=np.array([[0.,0.,0.],[1.,0.,0.],[0.,.5,0.]])
    xyz=np.stack([shape+[1.,1.,1.],shape+[7.5,7.5,7.5]])
    indicator=lambda x:int(len(eligible_acceptors(x,L,0,1,rlo,rhi,cmin))>0)
    current=indicator(xyz);saved=[];accepted=0
    for i in range(steps):
        trial,ratio,metadata=propose(xyz,L,0,rng,1,bias,rlo,rhi,cmin)
        bound=indicator(trial)
        logalpha=binding*(bound-current)+(ratio if corrected else 0.)
        if np.log(rng.random())<min(0.,logalpha):xyz=trial;current=bound;accepted+=1
        if i>=10000:saved.append(current)
    shell=4*np.pi*(rhi**3-rlo**3)/3;fraction=shell*(1-cmin)/2/L**3
    target=fraction*np.exp(binding)/(1-fraction+fraction*np.exp(binding))
    wrong_fraction=(1-bias)*fraction+bias
    wrong_target=wrong_fraction*np.exp(binding)/(1-wrong_fraction+wrong_fraction*np.exp(binding))
    measured=float(np.mean(saved));reference=target if corrected else wrong_target
    # Independent block uncertainty; this is not an estimate from HFIP.
    blocks=np.array(saved).reshape(90,1000).mean(axis=1)
    error=float(blocks.std(ddof=1)/np.sqrt(len(blocks)))
    qualified=bool(abs(measured-reference)<max(.0005,5*error))
    result=dict(binding_kBT=binding,Hastings_factor_applied=corrected,seed=seed,
                saved_steps=len(saved),geometric_association_probability=fraction,
                correct_target_probability=float(target),
                intentionally_wrong_target_if_factor_omitted=float(wrong_target),
                measured_association_probability=measured,block_standard_error=error,
                comparison_to_appropriate_analytic_distribution_passed=qualified,
                acceptance_fraction=accepted/steps)
    assert qualified,result
    return result


def main():
    root=HERE/'results/association-controls';assert not root.exists();root.mkdir(parents=True)
    start=time.monotonic();rows=[]
    for binding,corrected,seed in [(0.,True,20261030),(6.,True,20261031),(0.,False,20261032),(6.,False,20261033)]:
        row=run(binding,corrected,seed);rows.append(row)
        (root/'partial.json').write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
    result=dict(rows=rows,elapsed_seconds=time.monotonic()-start,
                accepted_angular_volume_sampler_control=True,
                reference_potential_not_HFIP=True,liquid_equilibrium_qualified=False,
                cleaning_solution_established=False)
    (root/'complete.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
