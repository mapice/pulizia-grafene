"""Scientific figure from completed derivative and same-start dynamics data."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
OUT=HERE.parent.parent/'outputs'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'cm',
                     'axes.unicode_minus':False,'font.size':10,'pdf.fonttype':42})
r=HERE/'results/virial-geometry-audit'
before=json.loads((r/'before.json').read_text());after=json.loads((r/'after.json').read_text())
b=HERE/'results/neural-bulk/HFIP_N16_seed20261016'
old=json.loads((b/'dynamics-pilot/NPT/observations.json').read_text())
new=json.loads((b/'pressure-corrected-NPT/observations.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(10.8,3.9))
for data,label,color in [(before,'Percorso precedente','#b96b43'),(after,'Derivata completa','#267d9b')]:
    rows=[x for x in data['checks'] if x['strain']=='isotropic']
    axes[0].loglog([x['step'] for x in rows],[abs(x['difference_eV_A3'])/6.241509074460763e-7 for x in rows],
                   '-o',label=label,color=color,ms=4)
axes[0].set(xlabel='Passo della deformazione della cella',ylabel='Scarto dalla derivata numerica (bar)',
            title='Controllo indipendente in doppia precisione')
axes[0].legend(frameon=False);axes[0].grid(alpha=.15)
for data,label,color,key in [(old,'Percorso precedente','#b96b43','conserved_extended_energy_eV'),
                             (new,'Derivata completa','#267d9b','extended_conserved_energy_eV')]:
    rows=[x for x in data if x['time_fs']<=50]
    e=np.array([x[key] for x in rows]);e-=e[0]
    axes[1].plot([x['time_fs'] for x in rows],e,'-o',ms=3,color=color,label=label)
axes[1].set(xlabel='Tempo (fs)',ylabel='Variazione dell’energia conservata estesa (eV)',
            title='Stessa partenza, controllo a pressione imposta')
axes[1].grid(alpha=.15);axes[1].legend(frameon=False)
fig.text(.5,.015,'HFIP, 16 molecole. Energia e forze preservate; controllo del metodo, senza qualificazione della pulizia o del liquido di equilibrio.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(OUT/'controllo-pressione-modello.pdf');fig.savefig(OUT/'controllo-pressione-modello.png',dpi=180)
print('Pressure/integration figure exported from completed evidence')
