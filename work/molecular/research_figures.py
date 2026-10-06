"""Standalone scientific figures from calculated data, without efficacy ranking."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
OUT=HERE.parent.parent/'outputs'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],
                     'mathtext.fontset':'cm','font.size':10,'axes.unicode_minus':False,
                     'pdf.fonttype':42,'ps.fonttype':42})

rotor=json.loads((HERE/'results/hfip-rotor/scan.json').read_text())['rows']
medium=json.loads((HERE/'results/polar-size-audit/M.json').read_text())['rotor']
gaff=json.loads((HERE/'results/hfip-rotor/gaff2.json').read_text())['rows']
fig,axes=plt.subplots(1,2,figsize=(11.8,4.2),gridspec_kw={'width_ratios':[1,1.25]})
ax=axes[0]
for rows,key,label,color in [(rotor,'total_relative_kJ_mol','PBE0-D3(BJ)/TZVP','#1d1d1d'),
                            (medium,'relative_energy_kJ_mol','MACE-POLAR-1-M','#1e6b8a'),
                            (gaff,'relative_energy_kJ_mol','GAFF2','#a96331')]:
    ax.plot([r['angle_degree'] for r in rows],[r[key] for r in rows],'-o',ms=3,lw=1.2,label=label,color=color)
ax.set(xlabel='Angolo di rotazione OH (gradi)',ylabel='Energia relativa (kJ/mol)',
       title='HFIP isolato: geometria fissata salvo H')
ax.set_xticks([-180,-90,0,90,180]);ax.legend(frameon=False,fontsize=9)
ax.grid(axis='y',alpha=.2)
labels=['HFIP','DMF','THF','Acido\nacetico','Anisolo']
values=[]
for name in ['HFIP','DMF','THF','AceticAcid','Anisole']:
    d=json.loads((HERE/'results/gaff-systems'/f'{name}_n4_seed20261006/interface-analysis.json').read_text())['statistics_second_half']
    values.append([d[key]['mean'] for key in ['backbone_and_side_methyl_G_LJ_kJ_mol',
                                           'ester_core_G_LJ_kJ_mol','ester_methyl_G_LJ_kJ_mol']])
values=np.array(values);ax=axes[1];bottom=np.zeros(5)
for j,(label,color) in enumerate([('Scheletro/metili','#787878'),('Nucleo estere','#bc7162'),('Metile estere','#e1be8b')]):
    ax.bar(np.arange(5),values[:,j],bottom=bottom,label=label,color=color,width=.65)
    bottom+=values[:,j]
ax.set_xticks(np.arange(5),labels);ax.set(ylabel='Interazione diretta con il grafene (kJ/mol)',
                    title='PMMA4: medie temporali di un breve campione')
ax.legend(frameon=False,fontsize=9,loc='lower left');ax.set_ylim(-100,5);ax.grid(axis='y',alpha=.2)
fig.text(.5,.018,'Modelli molecolari: nessuna energia libera di adsorbimento o efficacia di pulizia è dedotta da questi grafici.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.05,1,1])
fig.savefig(OUT/'controlli-molecolari-grafene.pdf')
fig.savefig(OUT/'controlli-molecolari-grafene.png',dpi=180)
plt.close(fig)
print('Standalone scientific figure saved')
