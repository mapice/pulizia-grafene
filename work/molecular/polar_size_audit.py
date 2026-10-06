"""Check whether the small-model field anomaly survives a larger model.

Neutral singlet inputs: code uses spin multiplicity1 (total_spin-1), not0.
All molecular geometries centred. Field-curvature sign does not depend on
whether the argument is the field or minus the field. No condensed-phase
or substrate qualification follows from successful molecular tests.
"""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[v] = '1'
os.nice(8)
from pathlib import Path
import argparse, gc, hashlib, json, time
import numpy as np
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
from ase.io import read
from mace.calculators import MACECalculator

HERE = Path(__file__).resolve().parent
HASHES = {'S': 'e4495612037b3b3312633182882a38a694ecac9ea0be2b9889ac0b2a84a99510',
          'M': 'fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--size', choices=HASHES, default='M')
    args = p.parse_args()
    path = HERE / 'sources' / f'MACE-POLAR-1-{args.size}.model'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == HASHES[args.size]
    torch.serialization.add_safe_globals([slice])
    model = torch.load(path, map_location='cpu', weights_only=False)
    calc = MACECalculator(models=model, device='cpu', default_dtype='float64',
                          model_type='PolarMACE', pbc_handling='realspace')
    del model
    gc.collect()
    out = HERE / 'results/polar-size-audit'
    out.mkdir(exist_ok=True)

    def evaluate(a, field):
        b = a.copy()
        b.info.update(charge=0, spin=1, external_field=list(field))
        b.calc = calc
        e = float(b.get_potential_energy())
        assert np.isfinite(e)
        return dict(energy_eV=e, field_argument_V_A=list(field),
                    density_dipole_eA=np.asarray(calc.results['dipole']).tolist(),
                    total_charge_e=float(np.sum(calc.results['charges'])))

    start = time.monotonic()
    checks = []
    for name in ['HFIP', 'DMF', 'THF']:
        a = read(HERE / 'results/clusters' / f'{name}.refined.xyz')
        a.positions -= a.positions.mean(0)
        zero = evaluate(a, [0., 0., 0.])
        axes = []
        for axis in range(3):
            delta = .01
            f = np.zeros(3)
            f[axis] = delta
            plus, minus = evaluate(a, f), evaluate(a, -f)
            alpha = -(plus['energy_eV'] + minus['energy_eV'] - 2 * zero['energy_eV']) / delta**2
            first = (plus['energy_eV'] - minus['energy_eV']) / (2 * delta)
            axes.append(dict(axis=axis, step_V_A=delta,
                             alpha_diagonal_eA2_V=alpha, energy_derivative_eA=first,
                             plus=plus, minus=minus))
        shifted = a.copy()
        shifted.positions += [3.5, -2., 4.]
        original_field = evaluate(a, [0., 0., .01])
        translated_field = evaluate(shifted, [0., 0., .01])
        error = translated_field['energy_eV'] - original_field['energy_eV']
        assert abs(error) < 1.e-7, error
        row = dict(molecule=name, zero=zero, axes=axes,
                   translation_energy_difference_eV=error,
                   all_diagonal_curvatures_positive=all(x['alpha_diagonal_eA2_V'] > 0 for x in axes))
        checks.append(row)
        print(name, [r['alpha_diagonal_eA2_V'] for r in axes], 'translation', error, flush=True)
    torsions = []
    for angle in range(-180, 180, 30):
        a = read(HERE / 'results/hfip-rotor' / f'angle{angle:+04d}.xyz')
        a.positions -= a.positions.mean(0)
        torsions.append(dict(angle_degree=angle, **evaluate(a, [0., 0., 0.])))
    emin = min(x['energy_eV'] for x in torsions)
    for r in torsions:
        r['relative_energy_kJ_mol'] = (r['energy_eV'] - emin) * 96.4853321233
    result = dict(model=f'MACE-POLAR-1-{args.size}', model_sha256=HASHES[args.size],
                  electronic_input='neutral charge, singlet multiplicity1; verified total_spin-1 in the installed official caller',
                  field_scope='Static electronic response at fixed nuclei, not bulk dielectric response',
                  fixed_gas_rotor_scope='Only the OH hydrogen is rotated; no solution torsional PMF',
                  response_checks=checks, rotor=torsions,
                  rotor_barrier_kJ_mol=max(x['relative_energy_kJ_mol'] for x in torsions),
                  elapsed_seconds=time.monotonic()-start,
                  condensed_phase_or_interface_qualification=False)
    (out / f'{args.size}.json').write_text(json.dumps(result, indent=2) + '\n')
    print('done', result['elapsed_seconds'], result['rotor_barrier_kJ_mol'], flush=True)


if __name__ == '__main__':
    main()
