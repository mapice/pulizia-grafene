"""Bounded sequential umbrella ladder, preserving pilots and both directions.

This produces biased molecular trajectories. A PMF is accepted only after
overlap, stationarity, autocorrelation and forward/reverse tests; the model
itself remains unqualified for a solvent efficacy ranking on PMMA950K.
"""
from pathlib import Path
import argparse, json, re, shutil, subprocess, time
import numpy as np

HERE = Path(__file__).resolve().parent
GMX = str(HERE / 'gromacs-env/bin/gmx')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--solvent', required=True, choices=['HFIP', 'DMF', 'Anisole'])
    p.add_argument('--direction', required=True, choices=['forward', 'reverse'])
    p.add_argument('--ps', type=int, default=500)
    p.add_argument('--k', type=float, default=750.)
    p.add_argument('--centers', default='.45,.55,.65,.75,.85,.95,1.05,1.15,1.25,1.35,1.45')
    p.add_argument('--max-minutes', type=float, default=240.)
    p.add_argument('--threads', type=int, choices=[1,2], default=2)
    p.add_argument('--pbc-previous-com',action='store_true')
    p.add_argument('--resume-failed-preparation',action='store_true')
    args = p.parse_args()
    base = HERE / 'results/gaff-systems' / f'{args.solvent}_n4_seed20261006'
    meta = json.loads((base / 'manifest.json').read_text())
    centers = sorted(float(x) for x in args.centers.split(','))
    assert len(set(centers)) == len(centers) and all(.35 <= x <= 1.5 for x in centers)
    if args.direction == 'reverse':
        centers.reverse()
        state = HERE / 'results/umbrella' / f'{args.solvent}_n4/z1.350_k2000'
        state_prefix = 'run'
        assert (state / 'complete.json').exists()
    else:
        state, state_prefix = base, 'production'
        assert json.loads((base / 'run-status.json').read_text())[-1]['exitcode'] == 0
    root = HERE / 'results/umbrella-ladders' / f'{args.solvent}_n4_k{args.k:g}_{args.ps}ps' / args.direction
    root.mkdir(parents=True, exist_ok=True)
    plan = dict(solvent=args.solvent, polymer_n=4, direction=args.direction,
                centers_nm=centers, force_constant_kJ_mol_nm2=args.k,
                ps_per_window=args.ps, initial_state=str(state),
                max_wall_minutes=args.max_minutes, threads=args.threads,
                scope='Conditional rigid-graphene oligomer model; not a cleaning prediction.',
                qualified_for_PMF=False, qualified_for_solvent_ranking=False)
    planfile = root / 'plan.json'
    if planfile.exists():
        old = json.loads(planfile.read_text())
        for key in ['solvent', 'direction', 'centers_nm', 'force_constant_kJ_mol_nm2', 'ps_per_window']:
            assert old[key] == plan[key], (key, old[key], plan[key])
    else:
        planfile.write_text(json.dumps(plan, indent=2) + '\n')
    rows, start_campaign = [], time.monotonic()
    for center in centers:
        w = root / f'z{center:.3f}'
        w.mkdir(exist_ok=True)
        if (w / 'complete.json').exists():
            rows.append(json.loads((w / 'complete.json').read_text()))
            state, state_prefix = w, 'run'
            continue
        if time.monotonic() - start_campaign > args.max_minutes * 60:
            (root / 'STOPPED-WALL-CAP.json').write_text(json.dumps(dict(completed=len(rows), cap_minutes=args.max_minutes)) + '\n')
            break
        if (w / 'run.tpr').exists():
            raise RuntimeError('Interrupted window preserved. Inspect and explicitly resume checkpoint: ' + str(w))
        if (w/'prepare.log').exists():
            assert args.resume_failed_preparation
            log=(w/'prepare.log').read_text()
            assert 'Set\n  the pull-pbc-ref-prev-step-com option to yes.' in log
            assert args.pbc_previous_com and not (w/'sample.log').exists()
            shutil.copy2(w/'prepare.log',w/'prepare-failed-pbc-reference.txt')
            shutil.copy2(w/'run.mdp',w/'input-before-pbc-tracking.mdp')
        for name in ['system.top', 'index.ndx']:
            shutil.copy2(base / name, w / name)
        mdp = re.sub(r'^nsteps\s*=\s*\d+\s*$', f'nsteps = {args.ps * 1000}',
                     (base / 'production.mdp').read_text(), flags=re.M)
        mdp += f'''
pull = yes
pull-ncoords = 1
pull-ngroups = 2
pull-group1-name = GRAPHENE
pull-group2-name = PMMA
pull-group1-pbcatom = 1
pull-group2-pbcatom = {meta['graphene_atoms'] + 14}
pull-coord1-type = umbrella
pull-coord1-geometry = direction-periodic
pull-coord1-groups = 1 2
pull-coord1-dim = N N Y
pull-coord1-vec = 0 0 1
pull-coord1-start = no
pull-coord1-init = {center:.9f}
pull-coord1-rate = 0
pull-coord1-k = {args.k:.9f}
pull-nstxout = 100
pull-nstfout = 100
'''
        if args.pbc_previous_com:
            mdp+='pull-pbc-ref-prev-step-com = yes\n'
        assert len(re.findall(r'^nsteps\s*=', mdp, re.M)) == 1
        (w / 'run.mdp').write_text(mdp)
        (w / 'start-state.json').write_text(json.dumps(dict(
            directory=str(state), prefix=state_prefix, center_nm=center,
            inherited_positions_and_velocities=True,
            pbc_reference_tracking_previous_COM=args.pbc_previous_com), indent=2) + '\n')
        commands = [('prepare', [GMX, 'grompp', '-f', 'run.mdp',
                                  '-c', str(state / (state_prefix + '.gro')),
                                  '-t', str(state / (state_prefix + '.cpt')),
                                  '-p', 'system.top', '-n', 'index.ndx', '-o', 'run.tpr']),
                    ('sample', [GMX, 'mdrun', '-deffnm', 'run', '-pf', 'pullf.xvg',
                                 '-px', 'pullx.xvg', '-ntmpi', '1', '-ntomp', str(args.threads),
                                 '-pin', 'off', '-nice', '8'])]
        started = time.monotonic()
        print(args.solvent, args.direction, center, 'started', flush=True)
        for stage, cmd in commands:
            with (w / (stage + '.log')).open('w') as f:
                r = subprocess.run([str(HERE / '.venv/bin/python'), str(HERE / 'guard_run.py'),
                                    '--mb', '1500', '--minutes', '35', '--'] + cmd,
                                   cwd=w, stdout=f, stderr=subprocess.STDOUT, timeout=2200)
            if r.returncode:
                raise RuntimeError('Window failed ' + str(w) + ' ' + stage)
        x = np.loadtxt(w / 'pullx.xvg', comments=['@', '#'])
        force = np.loadtxt(w / 'pullf.xvg', comments=['@', '#'])
        assert abs(x[-1, 0] - args.ps) < .11, x[-1, 0]
        discard = args.ps / 2
        tail, ftail = x[x[:, 0] >= discard], force[force[:, 0] >= discard]
        mean = float(tail[:, 1].mean())
        independently = args.k * (center - mean)
        measured = float(ftail[:, 1].mean())
        assert abs(measured - independently) < .01, (measured, independently)
        row = dict(solvent=args.solvent, direction=args.direction, center_nm=center,
                   k_kJ_mol_nm2=args.k, run_ps=args.ps, pilot_discard_ps=discard,
                   coordinate_mean_nm=mean, coordinate_sd_nm=float(tail[:, 1].std()),
                   spring_force_mean_kJ_mol_nm=measured,
                   independent_spring_force_kJ_mol_nm=independently,
                   elapsed_seconds=time.monotonic() - started,
                   pbc_reference_tracking_previous_COM=args.pbc_previous_com,
                   equilibrium_and_hidden_coordinate_convergence_not_assumed=True)
        (w / 'complete.json').write_text(json.dumps(row, indent=2) + '\n')
        rows.append(row)
        (root / 'windows.json').write_text(json.dumps(rows, indent=2) + '\n')
        print(row, flush=True)
        state, state_prefix = w, 'run'


if __name__ == '__main__':
    main()
