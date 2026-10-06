"""Compare actual opposite-history samples at identical bias, no PMF claim."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
root=HERE/'results/umbrella-ladders/HFIP_n4_k750_500ps'
out=HERE/'results/orientation-control'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'cm',
                     'axes.unicode_minus':False,'pdf.fonttype':42,'font.size':10})
keys=['CM_z_from_sheet_nm','Rg_normal_to_sheet_nm','longest_gyration_axis_normal_squared',
      'nearest_heavy_min_nm','smooth_heavy_contacts','polymer_graphene_LJ_energy_kJ_mol','HFIP_carbonyl_hbonds']
series={};summaries=[]
for direction in ['forward','reverse']:
    data=json.loads((root/direction/'z1.050/interface-analysis.json').read_text())
    tail=[r for r in data['rows'] if r['time_ps']>=250]
    series[direction]=np.array([[r[k] for k in keys] for r in tail])
    sub=[r for r in tail if 1.00<=r['CM_z_from_sheet_nm']<1.04]
    summaries.append(dict(direction=direction,frames=len(tail),sampling_interval_ps=1,
                          means={k:float(np.mean([r[k] for r in tail])) for k in keys},
                          matched_z_range_nm=[1.,1.04],matched_frames=len(sub),
                          matched_z_means={k:float(np.mean([r[k] for r in sub])) for k in keys} if sub else {},
                          frames_are_not_independent_samples=True))
z_edges=np.arange(.7,1.251,.01);r_edges=np.arange(.08,.341,.005)
hist=[np.histogram(series[d][:,0],z_edges)[0]/len(series[d]) for d in series]
joint=[np.histogram2d(series[d][:,0],series[d][:,1],bins=[z_edges,r_edges])[0]/len(series[d]) for d in series]
result=dict(bias_z_nm=1.05,bias_stiffness_kJ_mol_nm2=750,duration_ps_each=500,discarded_ps=250,
            histories=summaries,empirical_z_overlap_mass=float(np.minimum(*hist).sum()),
            empirical_z_rgz_overlap_mass=float(np.minimum(*joint).sum()),
            conclusion='Distance alone does not equilibrate the observed normal-orientation/contact mode in these samples',
            inference_scope='Single n4 rigid graphene/GAFF2 HFIP model; no true barrier, material cleaning or novelty established')
(out/'branch-comparison.json').write_text(json.dumps(result,indent=2)+'\n')
fig,axes=plt.subplots(1,3,figsize=(12,3.7))
labels={'forward':'Dal residuo aderente','reverse':'Dalla catena libera'}
colors={'forward':'#b56e3f','reverse':'#267d9b'}
for direction,x in series.items():
    axes[0].scatter(x[:,0],x[:,1],s=8,alpha=.35,label=labels[direction],color=colors[direction])
    axes[1].scatter(x[:,0],x[:,3],s=8,alpha=.35,color=colors[direction])
    axes[2].scatter(x[:,0],x[:,5],s=8,alpha=.35,color=colors[direction])
axes[0].set(ylabel='Estensione normale della catena (nm)',title='Orientazione diversa')
axes[1].set(ylabel='Distanza minima fra atomi pesanti e grafene (nm)',title='Contatto diverso')
axes[2].set(ylabel='Interazione diretta polimero–grafene (kJ/mol)',title='Energia diretta diversa')
for ax in axes:
    ax.set_xlabel('Altezza del centro di massa (nm)');ax.grid(alpha=.15)
axes[0].legend(frameon=False,fontsize=9)
fig.text(.5,.015,'Stesso vincolo sul centro di massa; ultimi 250 ps. Punti correlati: nessuna energia libera o efficacia di pulizia dedotta.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,1]);dest=HERE.parent.parent/'outputs'
fig.savefig(dest/'orientazione-e-contatti-grafene.pdf');fig.savefig(dest/'orientazione-e-contatti-grafene.png',dpi=180)
print(json.dumps(result,indent=2))
