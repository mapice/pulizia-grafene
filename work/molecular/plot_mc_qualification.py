"""Snapshot evidence: nonstationary HFIP chains and molecular internal energy."""
from pathlib import Path
import hashlib,json
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;OUT=HERE.parent.parent/'outputs'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'cm',
                    'axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42})
root=HERE/'results/neural-molecular-mc'
paths=[root/'HFIP_N16_seed20261026_1000moves/observations.json',
       root/'HFIP_N16_seed20261026_1000moves_continue2000/observations.json',
       root/'HFIP_N16_seed20261027_1000moves/observations.json']
data=[];provenance=[]
for path in paths:
    raw=path.read_bytes();rows=json.loads(raw);data.append(rows)
    provenance.append(dict(path=str(path.relative_to(HERE)),sha256=hashlib.sha256(raw).hexdigest(),
                           last_saved_attempt=rows[-1]['attempt'],terminal=(path.parent/'complete.json').exists()))
snap=HERE/'results/mc-internal-energy-audit/complete.json'
audit=json.loads(snap.read_text());molecules=audit['rows']
fig,axes=plt.subplots(1,2,figsize=(10.8,4.1))
dense=data[0]+data[1][1:]
for rows,label,color in [(dense,'Partenza a 1,607 g/cm³','#267d9b'),
                         (data[2],'Partenza a 1,000 g/cm³','#b96b43')]:
    axes[0].plot([r['attempt'] for r in rows],[r['instantaneous_density_g_cm3'] for r in rows],
                 label=label,color=color,lw=1.2)
axes[0].set(xlabel='Mosse tentate (nessun tempo fisico)',ylabel='Densità istantanea (g/cm³)',
            title='Due campionamenti ancora non confrontabili')
axes[0].legend(frameon=False,fontsize=9);axes[0].grid(alpha=.15)
selected=[molecules[i] for i in [0,2,3,4]]
values=[r['internal_distortion_kJ_mol_per_molecule'] for r in selected]
axes[1].bar(np.arange(4),values,color=['#267d9b','#267d9b','#b96b43','#b96b43'],alpha=.85,width=.65)
axes[1].set_xticks(np.arange(4),['Densa\niniziale','Densa\nsuccessiva','Diluita\niniziale','Diluita\nsuccessiva'])
axes[1].set(ylabel='Energia interna sopra il riferimento nel vuoto\n(kJ/mol per molecola)',
            title='Deformazioni interne ancora diverse')
for i,value in enumerate(values):axes[1].text(i,max(0,value)+.8,f'{value:.1f}',ha='center',fontsize=9)
axes[1].set_ylim(-3,43);axes[1].grid(axis='y',alpha=.15)
fig.text(.5,.014,'HFIP, 16 molecole, modello invariato. Configurazioni singole e serie transitorie: nessuna densità di equilibrio o efficacia di pulizia dedotta.',
         ha='center',fontsize=8.5)
fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(OUT/'qualificazione-hfip-liquido.pdf')
fig.savefig(OUT/'qualificazione-hfip-liquido.png',dpi=180)
record=dict(created_UTC=datetime.now(timezone.utc).isoformat(),observations=provenance,
            internal_audit_sha256=hashlib.sha256(snap.read_bytes()).hexdigest(),
            snapshots=selected,liquid_equilibrium_qualified=False,cleaning_solution_established=False,
            figure=str((OUT/'qualificazione-hfip-liquido.pdf').relative_to(HERE.parent.parent)))
(OUT/'qualificazione-hfip-liquido-dati.json').write_text(json.dumps(record,indent=2)+'\n')
print('Figure exported; terminal flags:',[(r['last_saved_attempt'],r['terminal']) for r in provenance])
