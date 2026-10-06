"""Numerical control of the full periodic molecular model on the Mac GPU.

Finite force/strain derivatives and two identical-start NVE trajectories
with step sizes 0.5 and 0.25 fs over 8 fs. This short control cannot give
density, dielectric constants, equilibrium structure or cleaning efficacy.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import json,time,hashlib
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float32);torch.mps.set_per_process_memory_fraction(.25)
from ase import units
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary
from ase.io import write
from mace.calculators import MACECalculator
from neural_structure import structure,SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=Path(__file__).resolve().parent


def main():
    root=HERE/'results/neural-liquid-numerics';root.mkdir(exist_ok=True)
    assert not (root/'summary.json').exists(),'Preserve completed experiment'
    parity=json.loads((HERE/'results/neural-liquid-bench/liquid64-parity.json').read_text())
    assert parity['force_relative_rms']<1e-4
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice])
    model=torch.load(path,map_location='cpu',weights_only=False).float()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='mps',default_dtype='float32')
    a=structure('liquid64');a.info.update(charge=0,spin=1,external_field=[0.,0.,0.]);a.calc=calc
    e=float(a.get_potential_energy());f=a.get_forces();stress=a.get_stress()
    checks=[]
    atom,axis=np.unravel_index(np.argmax(np.abs(f)),f.shape)
    delta=.002;energies=[]
    for sign in [-1,1]:
        b=a.copy();b.positions[atom,axis]+=sign*delta;b.calc=calc;energies.append(float(b.get_potential_energy()))
    numerical=-(energies[1]-energies[0])/(2*delta)
    checks.append(dict(type='atomic_force',atom=int(atom),axis=int(axis),step_A=delta,
                       analytical_eV_A=float(f[atom,axis]),numerical_eV_A=numerical,
                       relative_error=float(abs(numerical-f[atom,axis])/abs(f[atom,axis]))))
    molecular=f.reshape(64,12,3).sum(1);mol,axis=np.unravel_index(np.argmax(np.abs(molecular)),molecular.shape)
    delta=.01;energies=[]
    for sign in [-1,1]:
        b=a.copy();b.positions[mol*12:(mol+1)*12,axis]+=sign*delta;b.calc=calc;energies.append(float(b.get_potential_energy()))
    numerical=-(energies[1]-energies[0])/(2*delta)
    checks.append(dict(type='whole_molecule_force',molecule=int(mol),axis=int(axis),step_A=delta,
                       analytical_eV_A=float(molecular[mol,axis]),numerical_eV_A=numerical,
                       relative_error=float(abs(numerical-molecular[mol,axis])/abs(molecular[mol,axis]))))
    delta=.0001;energies=[]
    for sign in [-1,1]:
        b=a.copy();cell=a.cell.array.copy();cell[0,:]*=1+sign*delta;b.set_cell(cell,scale_atoms=True)
        b.calc=calc;energies.append(float(b.get_potential_energy()))
    numerical=(energies[1]-energies[0])/(2*delta*a.get_volume())
    checks.append(dict(type='xx_strain_derivative',strain_step=delta,analytical_eV_A3=float(stress[0]),
                       numerical_eV_A3=numerical,relative_error=float(abs(numerical-stress[0])/abs(stress[0]))))
    (root/'derivative-checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(checks,flush=True)
    assert all(c['relative_error']<.02 for c in checks),'Derivative discrepancy preserved; no dynamics accepted'
    a=structure('liquid64');a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    MaxwellBoltzmannDistribution(a,temperature_K=298.15,rng=np.random.default_rng(20261006));Stationary(a)
    initial_positions=a.positions.copy();initial_momenta=a.get_momenta().copy()
    write(root/'initial.xyz',a,format='extxyz')
    branches=[]
    for dt in [.5,.25]:
        b=a.copy();b.positions=initial_positions.copy();b.set_momenta(initial_momenta.copy());b.calc=calc
        rows=[];file=root/f'nve_dt{dt:.2f}fs.xyz';assert not file.exists()
        dyn=VelocityVerlet(b,dt*units.fs,logfile=None)
        start=time.monotonic()
        for step in range(round(8/dt)+1):
            ep=float(b.get_potential_energy());ek=float(b.get_kinetic_energy());et=ep+ek
            if not rows:e0=et
            row=dict(step=step,time_fs=step*dt,potential_eV=ep,kinetic_eV=ek,total_eV=et,
                     total_energy_change_eV=et-e0,temperature_K=float(b.get_temperature()),
                     net_momentum=b.get_momenta().sum(0).tolist())
            assert np.isfinite(b.positions).all() and np.isfinite(et)
            assert abs(et-e0)<.2,'Large numerical drift; preserve trajectory and stop'
            rows.append(row);write(file,b,format='extxyz',append=True)
            (root/f'nve_dt{dt:.2f}fs.json').write_text(json.dumps(rows,indent=2)+'\n')
            if step%4==0:print('NVE',dt,step,'dE',et-e0,flush=True)
            if step<round(8/dt):dyn.run(1)
        branches.append(dict(timestep_fs=dt,duration_fs=8,steps=round(8/dt),
                             max_absolute_total_energy_change_eV=float(max(abs(r['total_energy_change_eV']) for r in rows)),
                             final_total_energy_change_eV=rows[-1]['total_energy_change_eV'],
                             elapsed_seconds=time.monotonic()-start))
    result=dict(model_sha256=SHA,atoms=768,periodic_fixed_cell=True,initial_state='GAFF2-generated configuration, not M-model equilibrium',
                numerical_derivative_checks=checks,nve_branches=branches,
                MPS_memory_fraction_limit=.25,CPU_threads=1,scope='Numerical stability control only over 8fs',
                liquid_properties_or_interface_qualification=False,cleaning_solution_established=False)
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':main()
