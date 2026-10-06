"""Independent analytic controls for molecular-center NPT Monte Carlo.

Each group transformation x_i'=x_i+(s-1)COM has determinant s^3,
while preserving its internal coordinates. For symmetric proposals in
logV the acceptance exponent contains (N_molecules+1)*log(V'/V).
Rigid rotations and single OH rotations have unit Jacobian. These
controls concern the sampler, not HFIP or residue-removal predictions.
"""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(20261006)
    checks=[]
    for n in [3,7,12]:
        masses=rng.uniform(.8,20,n);w=masses/masses.sum()
        for ratio in [.8,1.,1.2]:
            s=ratio**(1/3)
            transform=np.eye(n)+(s-1)*np.outer(np.ones(n),w)
            A=np.kron(transform,np.eye(3))
            sign,logdet=np.linalg.slogdet(A)
            assert sign>0 and abs(logdet-np.log(ratio))<1e-12
            x=rng.normal(size=(n,3));c=np.average(x,axis=0,weights=masses)
            transformed=x+(s-1)*c
            assert np.max(np.abs((transformed-transformed[0])-(x-x[0])))<1e-12
            inverse=np.eye(n)+(1/s-1)*np.outer(np.ones(n),w)
            assert np.max(np.abs(inverse@transform-np.eye(n)))<1e-12
            checks.append(dict(atoms_in_group=n,volume_ratio=ratio,logdet_error=float(logdet-np.log(ratio))))
    # Independent known NPT distribution for ideal molecular centers:
    # beta*P=1, pi(V) proportional V^16 exp(-V): Gamma(shape17,scale1).
    n=16;logv=np.log(17.);values=[];accepted=0
    for step in range(220000):
        candidate=logv+rng.uniform(-.5,.5)
        exponent=-(np.exp(candidate)-np.exp(logv))+(n+1)*(candidate-logv)
        if np.log(rng.random())<min(0.,exponent):logv=candidate;accepted+=1
        if step>=20000 and step%20==0:values.append(np.exp(logv))
    values=np.array(values);mean=float(values.mean());variance=float(values.var())
    assert abs(mean-17)<.2 and abs(variance-17)<1.0,(mean,variance)
    result=dict(molecular_center_Jacobian_checks=checks,
                ideal_gas_NPT=dict(N_molecules=n,expected_mean_V=17.,expected_variance_V=17.,
                                  measured_mean_V=mean,measured_variance_V=variance,saved_samples=len(values),
                                  accept_fraction=accepted/220000),
                accepted_numerical_and_analytic_sampler_control=True,
                scope='Jacobian/ensemble control only, no molecular or cleaning result')
    root=HERE/'results/molecular-mc-controls';root.mkdir(exist_ok=True)
    (root/'control.json').write_text(json.dumps(result,indent=2)+'\n');print(result)


if __name__=='__main__':main()
