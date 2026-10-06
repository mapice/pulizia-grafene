"""Small OpenMM performance/sanity check; not a PMMA cleaning simulation."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
 os.environ[name]='2'
os.nice(8)
import openmm as mm
from openmm import unit
import numpy as np
import time,json,argparse,resource,sys
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--platform',choices=['CPU','OpenCL'],required=True)
args=p.parse_args()
system=mm.System()
nonbond=mm.NonbondedForce()
nonbond.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
nonbond.setCutoffDistance(1*unit.nanometer)
positions=np.array([(x,y,z) for x in range(6) for y in range(6) for z in range(6)],float)*.45
for _ in positions:
 system.addParticle(39.948)
 nonbond.addParticle(0,.34,.997)
system.addForce(nonbond)
system.setDefaultPeriodicBoxVectors(mm.Vec3(2.7,0,0),mm.Vec3(0,2.7,0),mm.Vec3(0,0,2.7))
integrator=mm.LangevinMiddleIntegrator(298.15*unit.kelvin,1/unit.picosecond,.001*unit.picoseconds)
integrator.setRandomNumberSeed(20261006)
platform=mm.Platform.getPlatformByName(args.platform)
props={'Threads':'2'} if args.platform=='CPU' else {'Precision':'mixed'}
start=time.monotonic()
context=mm.Context(system,integrator,platform,props)
context.setPositions(positions*unit.nanometer)
context.setVelocitiesToTemperature(298.15*unit.kelvin,20261006)
integrator.step(100)
warmup=time.monotonic()-start
start=time.monotonic()
integrator.step(2000)
seconds=time.monotonic()-start
state=context.getState(getEnergy=True,getPositions=True)
energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
assert np.isfinite(energy) and np.isfinite(xyz).all()
details={name:platform.getPropertyValue(context,name) for name in platform.getPropertyNames()}
peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)
report=dict(scope='engine benchmark, 216 generic Lennard-Jones particles, NOT sample/solvent prediction',platform=args.platform,platform_properties=details,atoms=216,steps=2000,dt_ps=.001,warmup_seconds=warmup,seconds=seconds,peak_rss_mb=peak,potential_kj_mol=energy)
(Path(__file__).resolve().parent/'results'/f'benchmark-{args.platform}.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
