"""Actual restrained separation of PMMA oligomer from rigid graphene in liquid.

Transparent harmonic-bias pilot with explicit contacts/hysteresis checks.
One-dimensional coordinate cannot certify convergence of hidden polymer
conformations. No950K/oxide-safety/cleaning claim.
"""
from pathlib import Path
import argparse,subprocess,time,json,shutil
import numpy as np
HERE=Path(__file__).resolve().parent
GMX=str(HERE/'gromacs-env/bin/gmx')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--solvent',default='HFIP')
    parser.add_argument('--centers',default='.60,.90,1.35')
    parser.add_argument('--ps',type=int,default=125);parser.add_argument('--k',type=float,default=2000)
    args=parser.parse_args();base=HERE/'results/gaff-systems'/f'{args.solvent}_n4_seed20261006'
    status=json.loads((base/'run-status.json').read_text())
    assert status[-1]['stage']=='production' and status[-1]['exitcode']==0
    root=HERE/'results/umbrella'/f'{args.solvent}_n4';root.mkdir(parents=True,exist_ok=True)
    meta=json.loads((base/'manifest.json').read_text())
    rows=[]
    for center in [float(x) for x in args.centers.split(',')]:
        folder=root/f'z{center:.3f}_k{args.k:g}';folder.mkdir(exist_ok=True)
        if (folder/'complete.json').exists():
            rows.append(json.loads((folder/'complete.json').read_text()));continue
        if (folder/'run.tpr').exists():raise RuntimeError('Existing interrupted window preserved; explicit checkpoint continuation required')
        for name in ['system.top','index.ndx']:shutil.copy2(base/name,folder/name)
        mdp=(base/'production.mdp').read_text().replace('nsteps = 100000',f'nsteps = {args.ps*1000}')
        mdp+=f'''
pull = yes
pull-ncoords = 1
pull-ngroups = 2
pull-group1-name = GRAPHENE
pull-group2-name = PMMA
pull-group1-pbcatom = 1
pull-group2-pbcatom = {meta['graphene_atoms']+14}
pull-coord1-type = umbrella
pull-coord1-geometry = direction-periodic
pull-coord1-groups = 1 2
pull-coord1-dim = N N Y
pull-coord1-vec = 0 0 1
pull-coord1-start = no
pull-coord1-init = {center:.9f}
pull-coord1-rate = 0.0
pull-coord1-k = {args.k:.9f}
pull-nstxout = 100
pull-nstfout = 100
'''
        (folder/'run.mdp').write_text(mdp)
        commands=[('prepare',[GMX,'grompp','-f','run.mdp','-c',str(base/'production.gro'),
                               '-t',str(base/'production.cpt'),'-p','system.top','-n','index.ndx',
                               '-o','run.tpr']),
                  ('sample',[GMX,'mdrun','-deffnm','run','-pf','pullf.xvg','-px','pullx.xvg',
                             '-ntmpi','1','-ntomp','2','-pin','off','-nice','8'])]
        start=time.monotonic()
        for stage,cmd in commands:
            with (folder/(stage+'.log')).open('w') as f:
                r=subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'guard_run.py'),
                                  '--mb','1500','--minutes','30','--']+cmd,
                                 cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=1900)
            if r.returncode:raise RuntimeError('Umbrella stage failed '+str(folder)+' '+stage)
        x=np.loadtxt(folder/'pullx.xvg',comments=['@','#'])
        f=np.loadtxt(folder/'pullf.xvg',comments=['@','#'])
        data=x[x[:,0]>=25];forces=f[f[:,0]>=25]
        record=dict(solvent=args.solvent,polymer_n=4,center_nm=center,k_kJ_mol_nm2=args.k,
                    run_ps=args.ps,discard_ps=25,
                    coordinate_mean_nm=float(data[:,1].mean()),coordinate_sd_nm=float(data[:,1].std()),
                    reported_pull_force_mean_kJ_mol_nm=float(forces[:,1].mean()),
                    independently_reconstructed_spring_force_kJ_mol_nm=args.k*(center-float(data[:,1].mean())),
                    elapsed_seconds=time.monotonic()-start,
                    scope='Single short biased window; not a converged free-energy profile or actual cleaning.')
        (folder/'complete.json').write_text(json.dumps(record,indent=2)+'\n')
        rows.append(record);(root/'windows.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(record,flush=True)

if __name__=='__main__':main()
