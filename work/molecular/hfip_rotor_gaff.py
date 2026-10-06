"""Independent OpenMM evaluation of the GAFF2 frozen gas OH rotation."""
from pathlib import Path
import json
import numpy as np
from ase.io import read
from openmm import app, unit, Context, VerletIntegrator, Platform

HERE = Path(__file__).resolve().parent
source = HERE / 'results/gaff2/HFIP'
top = app.AmberPrmtopFile(str(source / 'molecule.prmtop'))
coordinates = np.asarray(app.AmberInpcrdFile(str(source / 'molecule.inpcrd')).positions.value_in_unit(unit.angstrom))
z = np.array([a.element.atomic_number for a in top.topology.atoms()])
original = read(HERE / 'results/clusters/HFIP.refined.xyz')
assert np.array_equal(z, original.numbers)
err = float(abs(np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2) -
                np.linalg.norm(original.positions[:, None, :] - original.positions[None, :, :], axis=2)).max())
assert err < .005, err
system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
integrator = VerletIntegrator(.001 * unit.picoseconds)
context = Context(system, integrator, Platform.getPlatformByName('Reference'))
rows = []
for angle in range(-180, 180, 30):
    a = read(HERE / 'results/hfip-rotor' / f'angle{angle:+04d}.xyz')
    context.setPositions(a.positions * .1 * unit.nanometer)
    e = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    rows.append(dict(angle_degree=angle, energy_kJ_mol=e))
low = min(r['energy_kJ_mol'] for r in rows)
for r in rows:
    r['relative_energy_kJ_mol'] = r['energy_kJ_mol'] - low
result = dict(model='GAFF2/AM1-BCC neutral HFIP; original1-4LJ.5/Coulomb.833333333',
              scope='Frozen gas OH rotor; not solution conformer population, dielectric or cleaning prediction',
              atom_order_same_as_DFT=True, rounded_input_rigid_geometry_max_error_A=err,
              rows=rows, barrier_kJ_mol=max(r['relative_energy_kJ_mol'] for r in rows))
(HERE / 'results/hfip-rotor/gaff2.json').write_text(json.dumps(result, indent=2) + '\n')
print(result['barrier_kJ_mol'], [(x['angle_degree'], round(x['relative_energy_kJ_mol'], 3)) for x in rows])
