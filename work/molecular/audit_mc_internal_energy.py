"""Frozen snapshots: distinguish internal distortion from intermolecular energy.

This is a mechanical diagnostic, not a free-energy calculation. All fragment
energies use the same unchanged model, zero field, charge 0 and multiplicity 1.
Read complete XYZ frames once so concurrently running chains remain untouched.
"""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[v] = '1'
os.nice(8)
from pathlib import Path
from io import StringIO
import hashlib, json, time
import numpy as np
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float64)
from ase.io import read, write
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'results/mc-internal-energy-audit'
KJMOL = 96.4853321233


def snapshot(source, index, name):
    raw = source.read_bytes()
    lines = raw.decode().splitlines()
    n = int(lines[0]); assert n == 192
    complete = len(lines) // (n + 2)
    assert complete > 1
    actual = index if index >= 0 else complete + index
    block = '\n'.join(lines[actual*(n+2):(actual+1)*(n+2)]) + '\n'
    a = read(StringIO(block), format='extxyz')
    target = ROOT / (name + '.xyz')
    write(target, a, format='extxyz')
    return a, dict(source=str(source.relative_to(HERE)),
                   frame_index=actual, complete_frames_at_read=complete,
                   source_snapshot_sha256=hashlib.sha256(raw).hexdigest(),
                   saved_frame_sha256=hashlib.sha256(target.read_bytes()).hexdigest())


def main():
    assert not ROOT.exists(), 'Preserve already-completed or interrupted audit'
    ROOT.mkdir(parents=True)
    base = HERE / 'results/neural-molecular-mc'
    selected = [
        ('dense_start', base/'HFIP_N16_seed20261026_1000moves/frames.xyz', 0),
        ('dense_1000', base/'HFIP_N16_seed20261026_1000moves/frames.xyz', -1),
        ('dense_latest', base/'HFIP_N16_seed20261026_1000moves_continue2000/frames.xyz', -1),
        ('dilute_start', base/'HFIP_N16_seed20261027_1000moves/frames.xyz', 0),
        ('dilute_latest', base/'HFIP_N16_seed20261027_1000moves/frames.xyz', -1),
    ]
    frames = [(name, *snapshot(path, index, name)) for name, path, index in selected]
    path = HERE / 'sources/MACE-POLAR-1-M.model'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == SHA
    torch.serialization.add_safe_globals([slice])
    model = torch.load(path, map_location='cpu', weights_only=False).double()
    with torch.no_grad():
        model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model, 256); wrap_node_products(model, 16)
    calc = MACECalculator(models=model, model_type='PolarMACE',
                          device='cpu', default_dtype='float64')

    def energy(a):
        b = a.copy(); b.info.update(charge=0, spin=1, external_field=[0.,0.,0.])
        batch = calc._atoms_to_batch(b).to_dict()
        with torch.no_grad():
            result = model(batch, compute_force=False, compute_stress=False,
                           compute_virials=False, training=False)
        value = float(result['energy'].item()) * calc.energy_units_to_eV
        assert np.isfinite(value)
        return value

    mono = read(HERE/'results/neural-bulk/HFIP_N16_seed20261027_rho1.000/monomer.xyz')
    gas_reference = energy(mono)
    rows=[]; start=time.monotonic()
    for name, a, provenance in frames:
        length=float(a.cell[0,0]); xyz=a.positions.reshape(16,12,3)
        delta=xyz-xyz[:,:1,:]; delta-=length*np.rint(delta/length)
        assert np.linalg.norm(delta,axis=2).max()<length/2
        a.positions=(xyz[:,:1,:]+delta).reshape(-1,3)
        total=energy(a); fragments=[]
        for i in range(16):
            b=a[i*12:(i+1)*12]; b.pbc=False; b.set_cell(np.zeros((3,3)))
            b.positions-=b.positions.mean(axis=0)
            fragments.append(energy(b))
        a.info.update(charge=0,spin=1,external_field=[0.,0.,0.]); a.calc=calc
        force=a.get_forces(); stress=a.get_stress(voigt=False)
        # Configurational pressure only. The atom-based ideal kinetic term
        # does not make this unequilibrated snapshot an equation of state.
        from ase import units
        pconf=-float(np.trace(stress))/3/(1e5*units.Pascal)
        xyz=a.positions.reshape(16,12,3)
        oh=np.linalg.norm(xyz[:,10]-xyz[:,0],axis=1)
        co=np.linalg.norm(xyz[:,1]-xyz[:,0],axis=1)
        row=dict(name=name,provenance=provenance,
                 density_g_cm3=float(a.get_masses().sum()/(.602214076*a.get_volume())),
                 total_eV=total,isolated_fragment_eV=fragments,gas_reference_eV=gas_reference,
                 internal_distortion_kJ_mol_per_molecule=float((np.mean(fragments)-gas_reference)*KJMOL),
                 intermolecular_energy_kJ_mol_per_molecule=float((total-sum(fragments))/16*KJMOL),
                 total_above_gas_kJ_mol_per_molecule=float((total/16-gas_reference)*KJMOL),
                 OH_A_mean=float(oh.mean()),OH_A_range=[float(oh.min()),float(oh.max())],
                 CO_A_mean=float(co.mean()),CO_A_range=[float(co.min()),float(co.max())],
                 force_max_eV_A=float(np.linalg.norm(force,axis=1).max()),
                 configurational_pressure_bar=pconf,
                 elapsed_seconds=time.monotonic()-start)
        assert abs(row['total_above_gas_kJ_mol_per_molecule']-
                   row['internal_distortion_kJ_mol_per_molecule']-
                   row['intermolecular_energy_kJ_mol_per_molecule'])<1e-10
        rows.append(row)
        (ROOT/'partial.json').write_text(json.dumps(rows,indent=2)+'\n')
        print({k:row[k] for k in ['name','density_g_cm3','internal_distortion_kJ_mol_per_molecule',
                                  'intermolecular_energy_kJ_mol_per_molecule','configurational_pressure_bar']},flush=True)
    result=dict(model_sha256=SHA,device='CPU64',same_hamiltonian=True,
                atomic_reference_gauge_zero=True,rows=rows,
                scope='Frozen snapshot mechanical decomposition; no equilibrium chemical potential, entropy, EOS, solvent ranking or cleaning claim',
                liquid_equilibrium_qualified=False,cleaning_solution_established=False)
    (ROOT/'complete.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
