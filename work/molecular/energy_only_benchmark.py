"""Compare scalar-energy inference with official force/stress inference.

Same model and graph builder, compute_force=False only. Useful for exact
Metropolis moves; this benchmark by itself samples no liquid equilibrium.
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
from ase.io import read
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=Path(__file__).resolve().parent


def main():
    root=HERE/'results/energy-only-benchmark';root.mkdir(exist_ok=True)
    assert not (root/'comparison.json').exists()
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).float()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='mps',default_dtype='float32')
    source=HERE/'results/neural-bulk/HFIP_N16_seed20261016/dynamics-pilot/NVT/last-state.xyz'
    a=read(source);a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    rows=[]
    for perturbation in [0.,.01,-.01]:
        b=a.copy();b.positions[:12,0]+=perturbation;b.calc=calc
        start=time.monotonic();official=float(b.get_potential_energy());torch.mps.synchronize();time_full=time.monotonic()-start
        data=calc._atoms_to_batch(b).to_dict()
        start=time.monotonic()
        with torch.no_grad():
            output=model(data,compute_force=False,compute_stress=False,compute_virials=False,training=False)
            scalar=float(output['energy'].item())*calc.energy_units_to_eV
        torch.mps.synchronize();time_scalar=time.monotonic()-start
        error=scalar-official
        assert np.isfinite(scalar) and abs(error)<.0002
        rows.append(dict(translation_molecule0_x_A=perturbation,official_energy_eV=official,
                         scalar_energy_eV=scalar,difference_eV=error,
                         official_seconds=time_full,scalar_seconds=time_scalar))
        print(rows[-1],flush=True)
    result=dict(model_sha256=SHA,atoms=len(a),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                rows=rows,maximum_energy_difference_eV=float(max(abs(r['difference_eV']) for r in rows)),
                median_scalar_seconds=float(np.median([r['scalar_seconds'] for r in rows])),
                median_official_seconds=float(np.median([r['official_seconds'] for r in rows])),
                scope='Same Hamiltonian energy-path numerical control only; no equilibrium/ranking established')
    (root/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
