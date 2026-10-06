"""Analytic distribution and reversibility checks for random-length HMC.

Thirty independent harmonic modes stand in for classical intramolecular
coordinates. This diagnoses mixing of a sampler, not the HFIP potential.
"""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent


def leapfrog(q,p,omega,dt,steps):
    q=q.copy();p=p.copy();force=-(omega**2)*q
    for _ in range(steps):
        p+=.5*dt*force; q+=dt*p; force=-(omega**2)*q; p+=.5*dt*force
    return q,-p


def run(low,high,seed):
    rng=np.random.default_rng(seed); omega=np.linspace(.10,.70,30)
    dt=.25;q=np.zeros(30); values=[];accepted=0; early=[]
    for i in range(12000):
        p=rng.normal(size=30);old=.5*np.dot(omega*q,omega*q)+.5*np.dot(p,p)
        count=low if low==high else int(rng.integers(low,high+1))
        candidate,momentum=leapfrog(q,p,omega,dt,count)
        new=.5*np.dot(omega*candidate,omega*candidate)+.5*np.dot(momentum,momentum)
        if np.log(rng.random())<min(0.,old-new):q=candidate;accepted+=1
        potential=.5*np.dot(omega*q,omega*q)
        if i<50:early.append(potential)
        if i>=2000:values.append((omega*q)**2)
    values=np.array(values);mode_var=values.mean(axis=0)
    return dict(step_range_inclusive=[low,high],dt=.25,samples=len(values),
                expected_mode_variance=1.,mean_mode_variance=float(mode_var.mean()),
                maximum_mode_variance_error=float(np.max(np.abs(mode_var-1))),
                expected_total_potential_kBT=15.,
                measured_total_potential_kBT=float(.5*mode_var.sum()),
                first_50_total_potential_kBT=early,
                accepted_fraction=accepted/12000)


def main():
    rng=np.random.default_rng(20261028);omega=np.linspace(.10,.70,30)
    reversibility=[]
    for count in [4,12,19,32]:
        q=rng.normal(size=30);p=rng.normal(size=30)
        b,c=leapfrog(q,p,omega,.25,count);d,e=leapfrog(b,c,omega,.25,count)
        error=max(np.max(np.abs(d-q)),np.max(np.abs(e-p)))
        assert error<2e-12
        reversibility.append(dict(steps=count,phase_reversal_max_error=float(error)))
    rows=[run(4,4,20261028),run(12,32,20261029)]
    for row in rows:
        row['all_mode_precision_criterion_met']=(abs(row['mean_mode_variance']-1)<.04
                                                and row['maximum_mode_variance_error']<.13)
    result=dict(reversibility=reversibility,harmonic_modes=rows,
                accepted_random_length_sampler_control=True,
                scope='Analytic harmonic canonical distribution, not liquid/model physical qualification')
    root=HERE/'results/hmc-mixture-controls';root.mkdir(exist_ok=True)
    assert not (root/'control.json').exists()
    (root/'control.json').write_text(json.dumps(result,indent=2)+'\n')
    # Qualify the new mixture. Keep the short-kernel result as an explicitly
    # labelled mixing comparison, including failure to meet this precision.
    assert rows[1]['all_mode_precision_criterion_met'], rows[1]
    print({k:rows[1][k] for k in ['measured_total_potential_kBT','maximum_mode_variance_error','accepted_fraction']})


if __name__=='__main__':main()
