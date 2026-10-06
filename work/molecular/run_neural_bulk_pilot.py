"""Bounded actual neural-liquid pilot: fixed-volume then MTK variable-volume.

One small cell and 100fs per stage are preliminary checks, not equilibrium
qualification. The experimentally initialized density is never a result.
No latent-dipole dielectric inference is made.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float32);torch.mps.set_per_process_memory_fraction(.25)
from ase import units
from ase.io import read,write
from ase.optimize import FIRE
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary
from ase.md.nose_hoover_chain import NoseHooverChainNVT,IsotropicMTKNPT
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--molecules',type=int,choices=[16,32],default=16)
    p.add_argument('--seed',type=int,default=20261016);args=p.parse_args()
    base=HERE/'results/neural-bulk'/f'HFIP_N{args.molecules}_seed{args.seed}'
    prepared=json.loads((base/'prepared.json').read_text())
    assert hashlib.sha256((base/'initial.xyz').read_bytes()).hexdigest()==prepared['input_sha256']
    root=base/'dynamics-pilot';assert not root.exists(),'Preserve previous or interrupted experiment'
    root.mkdir()
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).float()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='mps',default_dtype='float32')
    a=read(base/'initial.xyz');a.info.update(charge=0,spin=1,external_field=[0.,0.,0.]);a.calc=calc
    plan=dict(model_sha256=SHA,molecules=args.molecules,atoms=len(a),timestep_fs=.25,
              fixed_volume_stage_fs=100,variable_volume_stage_fs=100,
              target_temperature_K=298.15,target_pressure_bar=1.,thermostat_damping_fs=25.,barostat_damping_fs=250.,
              NPT_integrator='ASE IsotropicMTKNPT; 1bar=1e5Pa, converted with ASE units.Pascal',
              MPS_memory_fraction_limit=.25,scope='Small-cell transient pilot only, not equilibrium properties or cleaning')
    (root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    optimizer=FIRE(a,logfile=str(root/'initial-minimization.txt'),maxstep=.05)
    relaxed=bool(optimizer.run(fmax=.15,steps=60))
    write(root/'relaxed-initial.xyz',a,format='extxyz')
    (root/'initial-minimization.json').write_text(json.dumps(dict(converged=relaxed,
              steps=int(optimizer.nsteps),fixed_initial_volume=True,equilibrium_not_implied=True),indent=2)+'\n')
    # Independent strain derivative before allowing a variable cell.
    stress=a.get_stress();delta=.0003;energies=[]
    for sign in [-1,1]:
        b=a.copy();cell=a.cell.array.copy();cell[0,:]*=1+sign*delta;b.set_cell(cell,scale_atoms=True);b.calc=calc
        energies.append(float(b.get_potential_energy()))
    numerical=(energies[1]-energies[0])/(2*delta*a.get_volume())
    error=abs(numerical-stress[0]);relative=error/max(abs(stress[0]),.0001)
    check=dict(analytical_xx_eV_A3=float(stress[0]),numerical_xx_eV_A3=float(numerical),
               relative_error_with_floor=float(relative),strain_step=delta)
    (root/'strain-check.json').write_text(json.dumps(check,indent=2)+'\n')
    assert relative<.02,'Virial check failed: variable-cell pilot rejected and retained'
    MaxwellBoltzmannDistribution(a,temperature_K=298.15,rng=np.random.default_rng(args.seed));Stationary(a)
    stages=[]
    for stage in ['NVT','NPT']:
        folder=root/stage;folder.mkdir()
        if stage=='NVT':
            dyn=NoseHooverChainNVT(a,timestep=.25*units.fs,temperature_K=298.15,tdamp=25*units.fs,logfile=None)
        else:
            dyn=IsotropicMTKNPT(a,timestep=.25*units.fs,temperature_K=298.15,pressure_au=1e5*units.Pascal,
                               tdamp=25*units.fs,pdamp=250*units.fs,logfile=None)
        rows=[];start=time.monotonic()
        for step in range(401):
            if step%20==0:
                length=float(min(a.cell.lengths()));temp=float(a.get_temperature())
                assert length>2*prepared['model_local_cutoff_A']+.2 and length<22
                assert 30<temp<1200 and np.isfinite(a.positions).all()
                e=float(a.get_potential_energy());fmax=float(np.linalg.norm(a.get_forces(),axis=1).max())
                assert np.isfinite(e) and fmax<50,'Numerical instability: retain pilot and stop'
                pressure=-float(np.trace(a.get_stress(voigt=False,include_ideal_gas=True)))/3
                row=dict(time_fs=step*.25,temperature_K=temp,volume_A3=float(a.get_volume()),
                         instantaneous_density_g_cm3=float(a.get_masses().sum()/(.602214076*a.get_volume())),
                         internal_pressure_bar=pressure/(1e5*units.Pascal),potential_eV=e,
                         max_force_eV_A=fmax,conserved_extended_energy_eV=float(dyn.get_conserved_energy()))
                rows.append(row);write(folder/'frames.xyz',a,format='extxyz',append=True)
                write(folder/'last-state.xyz',a,format='extxyz')
                (folder/'observations.json').write_text(json.dumps(rows,indent=2)+'\n')
                print(stage,row,flush=True)
            if step<400:dyn.run(1)
        result=dict(stage=stage,duration_fs=100,steps=400,elapsed_seconds=time.monotonic()-start,
                    initial_density_is_input=True,equilibrium_not_established=True,independent_replicas=1)
        (folder/'complete.json').write_text(json.dumps(result,indent=2)+'\n');stages.append(result)
    (root/'complete.json').write_text(json.dumps(dict(stages=stages,model_sha256=SHA,
              liquid_equilibrium_qualification=False,cleaning_solution_established=False),indent=2)+'\n')


if __name__=='__main__':main()
