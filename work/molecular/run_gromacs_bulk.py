"""Run one resource-bounded solvent qualification/control, never a lab action."""
from pathlib import Path
import subprocess,time,json,argparse

HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--scale14',type=float,default=1.0);p.add_argument('--seed',type=int,default=20261006);p.add_argument('--ps',type=int,default=500)
    args=p.parse_args()
    subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'build_gromacs_bulk.py'),'--scale14',str(args.scale14),'--seed',str(args.seed)],check=True)
    w=HERE/'results/gromacs'/f'bulk_scale{args.scale14:g}_seed{args.seed}'
    if (w/'npt.tpr').exists():raise RuntimeError('Existing run preserved; use explicit checkpoint continuation')
    mdp=w/'npt.mdp';mdp.write_text(mdp.read_text().replace('nsteps = 1000000',f'nsteps = {args.ps*1000}'))
    gmx=str(HERE/'gromacs-env/bin/gmx');rows=[]
    queue=[('min-prep',[gmx,'grompp','-f','min.mdp','-c','initial.gro','-p','system.top','-o','min.tpr']),('min',[gmx,'mdrun','-deffnm','min','-ntmpi','1','-ntomp','2','-pin','off','-nice','8']),('nvt-prep',[gmx,'grompp','-f','nvt.mdp','-c','min.gro','-p','system.top','-o','nvt.tpr']),('nvt',[gmx,'mdrun','-deffnm','nvt','-ntmpi','1','-ntomp','2','-pin','off','-nice','8']),('npt-prep',[gmx,'grompp','-f','npt.mdp','-c','nvt.gro','-t','nvt.cpt','-p','system.top','-o','npt.tpr']),('npt',[gmx,'mdrun','-deffnm','npt','-ntmpi','1','-ntomp','2','-pin','off','-nice','8'])]
    for stage,cmd in queue:
        start=time.monotonic();print('Starting',stage,flush=True)
        with (w/(stage+'.run.txt')).open('w') as f:
            r=subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'guard_run.py'),'--mb','1500','--minutes','30','--']+cmd,cwd=w,stdout=f,stderr=subprocess.STDOUT,timeout=1900)
        rows.append(dict(stage=stage,exitcode=r.returncode,elapsed_seconds=time.monotonic()-start))
        (w/'run-status.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows[-1],flush=True)
        if r.returncode:raise SystemExit(r.returncode)

if __name__=='__main__':main()
