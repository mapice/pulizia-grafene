"""Actual liquid HFIP internal-rotor distributions, with clear phase scope."""
from pathlib import Path
import subprocess,json,re
import numpy as np
from published_hfip import MOL,TYPES

HERE=Path(__file__).resolve().parent
GMX=str(HERE/'gromacs-env/bin/gmx')
nat=MOL.GetNumAtoms()
quartet=[TYPES.index(t)+1 for t in ['HO','OH','CT','HC']]

def xvg(p):
    return np.loadtxt(p,comments=['#','@'])

def main():
    summaries=[]
    for scale in [.5,1.]:
        p=HERE/'results/gromacs'/f'bulk_scale{scale:g}_seed20261006'
        index=p/'rotor.ndx'
        index.write_text('[ HFIP_OH_rotor ]\n'+'\n'.join(
            ' '.join(str(i+m*nat) for i in quartet) for m in range(125))+'\n')
        r=subprocess.run([GMX,'angle','-f','npt.xtc','-n','rotor.ndx',
                          '-type','dihedral','-all','-b','200','-dt','5',
                          '-od','rotor-hist.xvg','-ov','rotor-all.xvg'],
                         input='0\n',capture_output=True,text=True,cwd=p,timeout=60)
        (p/'rotor-analysis.txt').write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError(r.stderr[-1200:])
        values=xvg(p/'rotor-all.xvg')
        # GROMACS -all includes average followed by each individual torsion.
        if values.shape[1]!=127:raise RuntimeError('Unexpected angle columns '+str(values.shape))
        phi=values[:,2:].ravel()
        phi=(phi+180)%360-180
        circular=lambda target:np.abs((phi-target+180)%360-180)
        ap=float(np.mean(circular(180)<=30))
        syn=float(np.mean((circular(60)<=30)|(circular(-60)<=30)))
        hist,edges=np.histogram(phi,bins=np.arange(-180,185,5),density=True)
        row=dict(scale14=scale,first_ps=float(values[0,0]),last_ps=float(values[-1,0]),
                 frames=len(values),molecules=125,
                 observations=len(phi),observations_independent=False,
                 ap_fraction_within30deg=ap,
                 syn_fraction_within30deg=syn,
                 histogram_centers_degree=((edges[:-1]+edges[1:])/2).tolist(),
                 histogram_density_per_degree=hist.tolist(),
                 scope='Internal rotor in simulated liquid; no experimentally assigned population, no cleaning prediction',
                 conventional_quartet_one_based=quartet)
        (p/'rotor-summary.json').write_text(json.dumps(row,indent=2)+'\n')
        summaries.append(row);print(scale,'ap',ap,'syn',syn,'frames',len(values))
    (HERE/'results/hfip-liquid-rotors.json').write_text(json.dumps(summaries,indent=2)+'\n')

if __name__=='__main__':main()
