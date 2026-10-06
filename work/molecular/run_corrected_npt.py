"""Rerun the same NPT start after independently fixing cell derivatives.

Previous rejected variable-cell data are preserved. Start from the valid
NVT positions/momenta, with the same new MTK chain initialization as the
old NPT stage. 100fs is a transient/numerical check, not liquid validation.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import json,hashlib,time
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float32);torch.mps.set_per_process_memory_fraction(.25)
from ase import units
from ase.io import read,write
from ase.md.nose_hoover_chain import IsotropicMTKNPT
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products
HERE=Path(__file__).resolve().parent


def main():
    audit=json.loads((HERE/'results/virial-geometry-audit/comparison.json').read_text())
    assert audit['residual_after_small_step_bar']<.02 and audit['maximum_force_difference_eV_A']<1e-10
    base=HERE/'results/neural-bulk/HFIP_N16_seed20261016'
    source=base/'dynamics-pilot/NVT/last-state.xyz'
    root=base/'pressure-corrected-NPT';assert not root.exists()
    root.mkdir()
    a=read(source);assert a.has('momenta') and len(a)==192
    a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).float()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    a.calc=MACECalculator(models=model,model_type='PolarMACE',device='mps',default_dtype='float32')
    plan=dict(start_SHA256=hashlib.sha256(source.read_bytes()).hexdigest(),model_SHA256=SHA,
              same_positions_and_momenta_as_previous_NPT=True,cell_derivatives_corrected=True,
              dt_fs=.25,duration_fs=100,T_K=298.15,P_bar=1,tdamp_fs=25,pdamp_fs=250,
              liquid_equilibrium_not_established=True,MPS_memory_fraction_limit=.25)
    (root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    dyn=IsotropicMTKNPT(a,timestep=.25*units.fs,temperature_K=298.15,pressure_au=1e5*units.Pascal,
                       tdamp=25*units.fs,pdamp=250*units.fs,logfile=None)
    rows=[];start=time.monotonic()
    for step in range(401):
        if step%20==0:
            temp=float(a.get_temperature());length=float(min(a.cell.lengths()))
            assert 30<temp<1200 and 12.2<length<22 and np.isfinite(a.positions).all()
            e=float(a.get_potential_energy());stress=a.get_stress(voigt=False,include_ideal_gas=True)
            extended=float(dyn.get_conserved_energy())
            row=dict(time_fs=step*.25,temperature_K=temp,volume_A3=float(a.get_volume()),
                     instantaneous_density_g_cm3=float(a.get_masses().sum()/(.602214076*a.get_volume())),
                     internal_pressure_bar=-float(np.trace(stress))/(3e5*units.Pascal),
                     potential_eV=e,extended_conserved_energy_eV=extended,
                     extended_energy_change_eV=0. if not rows else extended-rows[0]['extended_conserved_energy_eV'])
            rows.append(row);write(root/'frames.xyz',a,format='extxyz',append=True)
            write(root/'last-state.xyz',a,format='extxyz')
            (root/'observations.json').write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
        if step<400:dyn.run(1)
    summary=dict(duration_fs=100,steps=400,elapsed_seconds=time.monotonic()-start,
                  maximum_extended_energy_change_eV=float(max(abs(r['extended_energy_change_eV']) for r in rows)),
                  initial_density_is_input_not_prediction=True,liquid_equilibrium_not_established=True,
                  cleaning_solution_established=False)
    (root/'complete.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
