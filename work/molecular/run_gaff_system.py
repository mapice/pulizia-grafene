"""Bounded actual dynamics for a prepared comparator system."""
from pathlib import Path
import argparse,subprocess,json,time
HERE=Path(__file__).resolve().parent
GMX=str(HERE/'gromacs-env/bin/gmx')

def main():
    p=argparse.ArgumentParser();p.add_argument('directory')
    p.add_argument('--from-stage',choices=['min-prep','nvt-prep','production-prep'],default='min-prep')
    p.add_argument('--threads',type=int,choices=[1,2],default=2)
    a=p.parse_args()
    w=Path(a.directory).resolve()
    if not w.is_relative_to(HERE/'results/gaff-systems'):raise RuntimeError('Wrong run directory')
    if (w/'production.tpr').exists():raise RuntimeError('Existing trajectory preserved; use explicit continuation')
    queue=[('min-prep',[GMX,'grompp','-f','min.mdp','-c','initial.gro','-p','system.top','-n','index.ndx','-o','min.tpr']),
           ('min',[GMX,'mdrun','-deffnm','min','-ntmpi','1','-ntomp','2','-pin','off','-nice','8']),
           ('nvt-prep',[GMX,'grompp','-f','nvt.mdp','-c','min.gro','-p','system.top','-n','index.ndx','-o','nvt.tpr']),
           ('nvt',[GMX,'mdrun','-deffnm','nvt','-ntmpi','1','-ntomp','2','-pin','off','-nice','8']),
           ('production-prep',[GMX,'grompp','-f','production.mdp','-c','nvt.gro','-t','nvt.cpt','-p','system.top','-n','index.ndx','-o','production.tpr']),
           ('production',[GMX,'mdrun','-deffnm','production','-ntmpi','1','-ntomp','2','-pin','off','-nice','8'])]
    if a.from_stage=='nvt-prep' and not (w/'min.gro').exists():raise RuntimeError('Missing minimized state')
    if a.from_stage=='production-prep' and not (w/'nvt.cpt').exists():raise RuntimeError('Missing equilibrated state')
    first=next(i for i,x in enumerate(queue) if x[0]==a.from_stage)
    rows=json.loads((w/'run-status.json').read_text()) if (w/'run-status.json').exists() else []
    for stage,cmd in queue[first:]:
        if '-ntomp' in cmd:cmd[cmd.index('-ntomp')+1]=str(a.threads)
        start=time.monotonic();print(w.name,'Starting',stage,flush=True)
        prior=w/(stage+'.run.txt')
        if prior.exists():
            k=1
            while (w/(stage+f'.previous{k}.txt')).exists():k+=1
            prior.rename(w/(stage+f'.previous{k}.txt'))
        with (w/(stage+'.run.txt')).open('w') as f:
            r=subprocess.run([str(HERE/'.venv/bin/python'),str(HERE/'guard_run.py'),
                              '--mb','1500','--minutes','30','--']+cmd,cwd=w,
                             stdout=f,stderr=subprocess.STDOUT,timeout=1900)
        rows.append(dict(stage=stage,exitcode=r.returncode,seconds=time.monotonic()-start))
        (w/'run-status.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(rows[-1],flush=True)
        if r.returncode:raise SystemExit(r.returncode)

if __name__=='__main__':main()
