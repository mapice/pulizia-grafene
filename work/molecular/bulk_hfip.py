"""Bulk HFIP qualification run with reproducible seed and resource bounds.

Published classical force field, not a graphene/residue simulation. Early
short runs are explicitly pilots, not converged dielectric measurements.
"""
import os
for var in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENMM_CPU_THREADS']:os.environ[var]='2'
os.nice(8)
import time,json,sys,resource,argparse
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from ase.io import read
import openmm as mm
from openmm import unit as u
from published_hfip import system,MOL

HERE=Path(__file__).resolve().parent
OUT=HERE/'results/bulk'
OUT.mkdir(parents=True,exist_ok=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--scale14',type=float,default=.5);p.add_argument('--ps',type=float,default=100);p.add_argument('--seed',type=int,default=20261006);p.add_argument('--restart',default=None);p.add_argument('--side',type=int,default=4)
    args=p.parse_args();n=args.side**3;nat=MOL.GetNumAtoms();temp=298.15
    rng=np.random.default_rng(args.seed)
    mass=168.0404
    # NA * 1e-21 cm^3/nm^3 = 602.214076, for g/mol and g/cm^3.
    box=(n*mass/(602.214076*.8))**(1/3)
    s,charges,mass=system(n,box,args.scale14,cutoff=.9 if args.side==4 else 1.0)
    barostat=mm.MonteCarloBarostat(1*u.bar,temp*u.kelvin,25);barostat.setRandomNumberSeed(args.seed+1);s.addForce(barostat)
    integrator=mm.LangevinMiddleIntegrator(temp*u.kelvin,1/u.picosecond,.001*u.picosecond);integrator.setRandomNumberSeed(args.seed);integrator.setConstraintTolerance(1e-6)
    context=mm.Context(s,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'2'})
    tag=f'hfip_scale{args.scale14:g}_seed{args.seed}'
    if args.restart:
        context.loadCheckpoint(Path(args.restart).read_bytes())
    else:
        a=read(HERE/'results/clusters/HFIP.refined.xyz');local=a.positions*.1;local-=local.mean(axis=0)
        positions=[]
        for i in range(args.side):
            for j in range(args.side):
                for k in range(args.side):positions.extend(Rotation.random(random_state=rng).apply(local)+(np.array([i,j,k])+.5)*box/args.side)
        context.setPositions(np.array(positions)*u.nanometer)
        print('Context initialized, minimizing',flush=True)
        mm.LocalEnergyMinimizer.minimize(context,50*u.kilojoule_per_mole/u.nanometer,100)
        print('Minimization finished',flush=True)
        context.setVelocitiesToTemperature(temp*u.kelvin,args.seed)
    rows=[];start=time.monotonic();steps=int(round(args.ps/.001));chunk=1000
    for done in range(0,steps,chunk):
        integrator.step(min(chunk,steps-done))
        st=context.getState(getEnergy=True,getPositions=True,enforcePeriodicBox=True)
        vec=np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer));v=float(np.linalg.det(vec))
        rho=n*mass/(602.214076*v)
        xyz=st.getPositions(asNumpy=True).value_in_unit(u.nanometer).reshape(n,nat,3)
        # Minimum-image molecule unwrapping before constructing the total dipole.
        anchor=xyz[:,:1,:];frac=(xyz-anchor)@np.linalg.inv(vec);xyz=anchor+(frac-np.rint(frac))@vec
        dip=np.sum(xyz*charges[None,:,None],axis=(0,1)) # e nm
        e=float(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole))
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)
        if not np.isfinite(e) or not np.isfinite(xyz).all():raise RuntimeError('Nonfinite MD state')
        if rss>3000:raise MemoryError('3 GB MD RSS cap')
        if time.monotonic()-start>1800:raise TimeoutError('30 min pilot cap; checkpoint saved after preceding chunk')
        row=dict(time_ps=(done+min(chunk,steps-done))*.001,density_g_cm3=rho,volume_nm3=v,energy_kJ_mol=e,dipole_e_nm=dip.tolist(),peak_rss_MB=rss,elapsed_seconds=time.monotonic()-start)
        rows.append(row)
        if len(rows)%10==0:
            (OUT/(tag+'.checkpoint')).write_bytes(context.createCheckpoint())
            (OUT/(tag+'.json')).write_text(json.dumps(dict(scope='Bulk HFIP qualification pilot; not cleaning efficacy; dielectric not converged in short trajectory',scale14=args.scale14,molecules=n,temperature_K=temp,seed=args.seed,source_doi='10.1002/cphc.202100620',rows=rows),indent=2)+'\n')
            print(json.dumps(row),flush=True)
    (OUT/(tag+'.checkpoint')).write_bytes(context.createCheckpoint())
    np.savez_compressed(OUT/(tag+'.final.npz'),positions=st.getPositions(asNumpy=True).value_in_unit(u.nanometer),box=vec)
    print('Finished',tag,rows[-1],flush=True)

if __name__=='__main__':main()
