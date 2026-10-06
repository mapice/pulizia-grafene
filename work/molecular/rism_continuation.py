"""Documented density/charge continuation for stiff organic DRISM equations.

Amber26 manual8.3.4.3 permits using easier artificial liquids to seed the
target liquid. Only the last, exact target state can provide an accepted
susceptibility. Bootstrap states are never cleaning/solvation predictions.
"""
from pathlib import Path
import argparse, hashlib, json, re, shutil, subprocess, time
import numpy as np

HERE = Path(__file__).resolve().parent


def scale_charges(text, factor):
    before, rest = text.split('%FLAG CHG\n', 1)
    block, after = rest.split('%FLAG LJEPSILON\n', 1)
    lines = block.splitlines()
    numbers = [float(x) * factor for x in ' '.join(lines[1:]).split()]
    formatted = '\n'.join(''.join(f'{x:16.8E}' for x in numbers[j:j+5])
                           for j in range(0, len(numbers), 5))
    return before + '%FLAG CHG\n' + lines[0] + '\n' + formatted + '\n%FLAG LJEPSILON\n' + after


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--solvent', choices=['HFIP', 'DMF'], default='HFIP')
    p.add_argument('--max-minutes', type=float, default=20.)
    args = p.parse_args()
    src = HERE / 'results/rism' / args.solvent
    meta = json.loads((src / 'manifest.json').read_text())
    target_mdl, target_input = (src / 'solvent.mdl').read_text(), (src / 'liquid.inp').read_text()
    root = HERE / 'results/rism' / (args.solvent + '-continuation')
    root.mkdir(exist_ok=True)
    if (root / 'attempts.json').exists():
        raise RuntimeError('Previous branch preserved; explicit review/restart required')
    # These are numerical seeds, without a claim to an experimentally
    # realisable dielectric response. The final state restores every input.
    states = [(.1, r) for r in [.1, .25, .5, .75, 1.]]
    states += [(q, 1.) for q in [.2, .3, .4, .5, .6, .7, .8, .9, .95, 1.]]
    attempts, previous, started = [], None, time.monotonic()
    for i, (qscale, rho_fraction) in enumerate(states):
        if time.monotonic() - started > 60 * args.max_minutes:
            break
        out = root / f'{i:02d}_q{qscale:.2f}_rho{rho_fraction:.2f}'
        out.mkdir(exist_ok=True)
        if (out / 'solve.log').exists():
            raise RuntimeError('Existing trial preserved')
        final = (qscale == 1. and rho_fraction == 1.)
        (out / 'solvent.mdl').write_text(target_mdl if final else scale_charges(target_mdl, qscale))
        input_text = target_input
        input_text = input_text.replace('MDIIS_DEL=0.3', 'MDIIS_DEL=0.1, MDIIS_RESTART=1000')
        input_text = input_text.replace('MAXSTEP=10000', 'MAXSTEP=1500').replace('KSAVE=-1', 'KSAVE=100')
        dielectric = meta['dielectric_imposed_from_experiment'] if final else (
            1 + (meta['dielectric_imposed_from_experiment']-1) * qscale**2 * rho_fraction)
        input_text = re.sub(r'DIEPS=[\d.]+', f'DIEPS={dielectric:.12f}', input_text)
        density = meta['molarity_mol_L'] * rho_fraction
        input_text = re.sub(r'DENSITY=[\d.]+', f'DENSITY={density:.12f}', input_text)
        input_text = input_text.replace('TOLERANCE=1.e-10', 'TOLERANCE=1.e-10' if final else 'TOLERANCE=1.e-8')
        (out / 'liquid.inp').write_text(input_text)
        if previous is not None:
            shutil.copy2(previous / 'liquid.sav', out / 'liquid.sav')
        with (out / 'solve.log').open('w') as f:
            r = subprocess.run([str(HERE / '.venv/bin/python'), str(HERE / 'guard_run.py'),
                                '--mb', '300', '--minutes', '2', '--', '/usr/bin/nice', '-n', '8',
                                str(HERE / 'amber-env/bin/rism1d'), 'liquid'],
                               cwd=out, stdout=f, stderr=subprocess.STDOUT, timeout=140)
        log = (out / 'solve.log').read_text()
        matches = re.findall(r'Res=\s*([\d.E+\-]+)', log)
        residual = float(matches[-1]) if matches else None
        row = dict(stage=i, charge_scale=qscale, density_fraction=rho_fraction,
                   dielectric_input=dielectric, exitcode=r.returncode, residual=residual,
                   xvv_written=(out / 'liquid.xvv').exists(),
                   exact_target_state=final, earlier_converged_restart=previous is not None)
        attempts.append(row)
        (root / 'attempts.json').write_text(json.dumps(attempts, indent=2) + '\n')
        print(row, flush=True)
        if r.returncode or not (out / 'liquid.xvv').exists() or not (out / 'liquid.sav').exists():
            (out / 'NOT-CONVERGED.json').write_text(json.dumps(row, indent=2) + '\n')
            break
        previous = out
        if final:
            assert hashlib.sha256((out / 'solvent.mdl').read_bytes()).hexdigest() == hashlib.sha256((src / 'solvent.mdl').read_bytes()).hexdigest()
    accepted = bool(attempts and attempts[-1]['exact_target_state'] and
                    attempts[-1]['exitcode'] == 0 and attempts[-1]['xvv_written'])
    (root / 'status.json').write_text(json.dumps(dict(
        target_solver_converged=accepted, source_parameters=meta,
        previous_unconverged_trials_not_used=True,
        target_numerical_and_physical_qualification_still_required=True,
        primary_manual='Amber26 sections8.3.4.2-3',
        elapsed_seconds=time.monotonic()-started), indent=2) + '\n')


if __name__ == '__main__':
    main()
