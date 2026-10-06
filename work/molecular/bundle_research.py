"""Snapshot completed research only; never bundle live trajectories or runtimes."""
from pathlib import Path
import datetime,hashlib,json,zipfile
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
selected={}

def add(path,label=None):
    if path.is_file():
        selected[label or 'molecular/'+str(path.relative_to(HERE))]=path

for p in HERE.glob('*.py'):add(p)
for name in ['LEGGIMI.txt','requirements-lock.txt','gromacs-explicit-lock.txt','amber-explicit-lock.txt']:
    add(HERE/name)
for name in ['model_cleaning.py','model_review.py','verify_models.py','terzo_ancoraggi.py','atomistic-campaign.json']:
    p=PROJECT/'work'/name
    if p.exists():selected['theory/'+name]=p

skip_dirs={'__pycache__','invalid-volume-factor1000','HFIP-original-nonneutral'}
allowed={'.json','.txt','.log','.md','.xyz','.mol2','.sdf','.pdb','.prmtop','.inpcrd',
         '.itp','.top','.atp','.gro','.g96','.mdp','.ndx','.tpr','.cpt','.edr','.xtc',
         '.xvg','.npz','.sav','.mdl','.inp','.xvv','.therm','.self.test','.dat'}
root=HERE/'results'
for p in root.rglob('*'):
    if not p.is_file() or any(part in skip_dirs for part in p.relative_to(root).parts):continue
    if any('invalid-input' in part or 'failed-pack' in part for part in p.relative_to(root).parts):continue
    if p.name.startswith(('analysis-whole','input-parsed')):continue
    if p.name.endswith('.scf.log') and p.parent.name.startswith('dft-matched'):
        # Only accepted SCF logs, including complete constituents before a
        # five-point comparison closes. Do not include the aborted dense run.
        record=p.with_name(p.name.replace('.scf.log','.json'))
        if not record.exists() or not json.loads(record.read_text()).get('SCF_converged',False):continue
    parts=p.relative_to(root).parts
    if parts[0].startswith('cf4-longrange-audit') and not (root/parts[0]/'complete.json').exists():continue
    if parts[0]=='gaff-systems':
        run_dir=root/parts[0]/parts[1]
        status=run_dir/'run-status.json'
        if not status.exists():continue
        rows=json.loads(status.read_text())
        if not rows or rows[-1]['stage']!='production' or rows[-1]['exitcode']!=0:continue
    if parts[0] in ['umbrella','umbrella-ladders']:
        window=next((q for q in p.parents if q.name.startswith('z0.') or q.name.startswith('z1.')),None)
        if window is not None and not (window/'complete.json').exists():continue
    if parts[0]=='orientation-pilots':
        branch=next((q for q in p.parents if q.name in ['forward','reverse']),None)
        if branch is not None and not (branch/'complete.json').exists():continue
    if parts[0]=='neural-liquid-numerics' and not (root/parts[0]/'summary.json').exists():continue
    if parts[0]=='neural-bulk' and 'dynamics-pilot' in parts:
        pilot=next(q for q in p.parents if q.name=='dynamics-pilot')
        if not (pilot/'complete.json').exists() and not (pilot/'rejected-pressure-audit.json').exists():continue
    if parts[0]=='neural-bulk' and 'pressure-corrected-NPT' in parts:
        branch=next(q for q in p.parents if q.name=='pressure-corrected-NPT')
        if not (branch/'complete.json').exists():continue
    if parts[0]=='neural-molecular-mc':
        branch=root/parts[0]/parts[1]
        if not (branch/'complete.json').exists():
            interruption=branch/'interrupted.json'
            if not interruption.exists() or not json.loads(interruption.read_text()).get('authoritatively_terminal',False):continue
    if p.suffix in allowed and p.stat().st_size<=50_000_000:add(p)

# Source metadata only: not full copyrighted articles or model weights.
for p in (HERE/'sources').rglob('*'):
    if p.is_file() and p.suffix in ['.json','.yaml'] and p.stat().st_size<1_000_000:add(p)
for p in (HERE/'patches').glob('*.json'):add(p)
for p in (HERE/'logs').glob('*.txt'):
    if p.stat().st_size<2_000_000:add(p)

# Distribute only licensed vendor code needed for a repeatable parameter
# inference. No install scripts, datasets or executables are included.
vendor=HERE/'vendor/byteff-pol'
for folder in ['byteff2','submodules/bytemol/bytemol']:
    for p in (vendor/folder).rglob('*'):
        if p.is_file() and p.suffix in ['.py','.csv','.json','.yaml'] and '__pycache__' not in p.parts:
            add(p)
for p in vendor.glob('*LICENSE*'):add(p)

for name in ['pulizia-grafene-pmma-ppc.tex','pulizia-grafene-pmma-ppc.pdf','controlli-molecolari-grafene.pdf','orientazione-e-contatti-grafene.pdf','controllo-pressione-modello.pdf','qualificazione-hfip-liquido.pdf','qualificazione-hfip-liquido-dati.json']:
    selected['paper/'+name]=PROJECT/'outputs'/name
entries=[]
manifest=dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),revision=7,
              terminal_data_snapshot=True,confirmed_interrupted_data_explicitly_labelled=True,live_trajectories_excluded=True,
              model_weights_and_runtimes_excluded=True,experimental_solution_validated=False,
              files=entries)
target=PROJECT/'work/grafene-modelli-revisione-7-new.zip'
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name,p in sorted(selected.items()):
        # Read once: a log/status append cannot desynchronize its digest.
        data=p.read_bytes()
        entries.append(dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
        z.writestr(name,data)
    z.writestr('MANIFEST.json',json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
with zipfile.ZipFile(target) as z:
    assert z.testzip() is None
    stored=json.loads(z.read('MANIFEST.json'))
    assert len(stored['files'])==len(selected)
    for row in stored['files']:
        assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256'],row['path']
target.replace(PROJECT/'outputs/grafene-modelli-e-sorgenti.zip')
target=PROJECT/'outputs/grafene-modelli-e-sorgenti.zip'
print(json.dumps(dict(archive=str(target),files=len(selected),bytes=target.stat().st_size,
                      sha256=hashlib.sha256(target.read_bytes()).hexdigest(),integrity_verified=True)))
