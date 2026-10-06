"""Actual two-coordinate pilot, started independently from the two branches.

Only tests access to the newly identified orientation coordinate; no PMF,
cleaning ranking, long-chain result or real-substrate safety is inferred.
"""
from pathlib import Path
import argparse,json,re,shutil,subprocess,time,os
import numpy as np
from orientation_cv import specification

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--ps',type=int,default=100)
    p.add_argument('--z',type=float,default=1.05);p.add_argument('--rgz',type=float,default=.22)
    args=p.parse_args();assert 50<=args.ps<=500 and .6<=args.z<=1.3 and .12<=args.rgz<=.35
    check=json.loads((HERE/'results/orientation-control/validation.json').read_text())
    assert check['coordinate_and_first_derivative_validation_passed']
    base=HERE/'results/gaff-systems/HFIP_n4_seed20261006'
    root=HERE/'results/orientation-pilots'/f'HFIP_z{args.z:.3f}_rgz{args.rgz:.3f}_{args.ps}ps'
    root.mkdir(parents=True,exist_ok=True)
    plan=dict(z_center_nm=args.z,rg_normal_center_nm=args.rgz,stiffnesses_kJ_mol_nm2=[750,2000],
              duration_ps=args.ps,threads=1,scope='Implementation/sampling pilot in unqualified n4 rigid-sheet model',
              desorption_equilibrium_or_cleaning_established=False)
    (root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    env=os.environ.copy();env['OMP_NUM_THREADS']='1'
    env['PLUMED_KERNEL']=str(HERE/'gromacs-env/lib/libplumedKernel.dylib')
    for direction in ['forward','reverse']:
        folder=root/direction;folder.mkdir(exist_ok=True)
        if (folder/'complete.json').exists():continue
        if (folder/'run.tpr').exists():raise RuntimeError('Interrupted pilot preserved, explicit checkpoint resume required')
        state=HERE/'results/umbrella-ladders/HFIP_n4_k750_500ps'/direction/'z1.050'
        assert (state/'complete.json').exists()
        for f in ['system.top','index.ndx']:shutil.copy2(base/f,folder/f)
        mdp=re.sub(r'^nsteps\s*=\s*\d+\s*$',f'nsteps = {1000*args.ps}',(base/'production.mdp').read_text(),flags=re.M)
        assert 'pull =' not in mdp
        (folder/'run.mdp').write_text(mdp)
        (folder/'plumed.dat').write_text(specification(file='CV',stride=100,bias=(args.z,args.rgz,750,2000)))
        (folder/'start-state.json').write_text(json.dumps(dict(directory=str(state),positions_and_velocities_inherited=True),indent=2)+'\n')
        commands=[('prepare',[str(HERE/'gromacs-env/bin/gmx'),'grompp','-f','run.mdp','-c',str(state/'run.gro'),
                              '-t',str(state/'run.cpt'),'-p','system.top','-n','index.ndx','-o','run.tpr']),
                  ('sample',[str(HERE/'gromacs-env/bin/gmx'),'mdrun','-deffnm','run','-plumed','plumed.dat',
                             '-ntmpi','1','-ntomp','1','-pin','off','-nice','8'])]
        start=time.monotonic();print('Started orientation pilot',direction,flush=True)
        for label,cmd in commands:
            with (folder/(label+'.log')).open('w') as log:
                r=subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'guard_run.py'),'--mb','1500','--minutes','20','--']+cmd,
                                 cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1250)
            if r.returncode:raise RuntimeError(f'Pilot {direction} {label} failed; logs retained')
        values=np.loadtxt(folder/'CV',comments='#')
        assert abs(values[-1,0]-args.ps)<.11
        tail=values[values[:,0]>=args.ps/2]
        result=dict(direction=direction,run_ps=args.ps,pilot_discard_ps=args.ps/2,
                    mean_z_nm=float(tail[:,1].mean()),mean_rgz_nm=float(tail[:,2].mean()),
                    elapsed_seconds=time.monotonic()-start,conditional_single_start_sampling_only=True)
        (folder/'complete.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':main()
