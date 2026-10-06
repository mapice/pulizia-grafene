"""Independent thermodynamic audit of generated GAFF2 solvent trajectories."""
from pathlib import Path
import subprocess,argparse,json,re
import numpy as np
from analyze_bulk import xvg
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('--start-ps',type=float,default=100);a=p.parse_args()
    w=Path(a.directory).resolve();meta=json.loads((w/'manifest.json').read_text())
    assert meta['polymer_n']==0
    g=str(HERE/'gromacs-env/bin/gmx')
    r=subprocess.run([g,'energy','-f','production.edr','-o','thermo.xvg'],
                     input='Temperature\nDensity\nVolume\n0\n',cwd=w,capture_output=True,text=True,timeout=30)
    (w/'thermo-analysis.txt').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Thermo analysis failed')
    r=subprocess.run([g,'dipoles','-s','production.tpr','-f','production.xtc',
                      '-b',str(a.start_ps),'-o','dipoles.xvg','-eps','epsilon.xvg',
                      '-temp','298.15','-nopairs'],input='0\n',cwd=w,capture_output=True,text=True,timeout=30)
    (w/'dipoles-analysis.txt').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Dipole analysis failed')
    values,columns=xvg(w/'thermo.xvg');tail=values[values[:,0]>=a.start_ps]
    dip,_=xvg(w/'dipoles.xvg');m=dip[:,1:4]
    var=float((m*m).sum(axis=1).mean()-(m.mean(axis=0)**2).sum())
    vol=float(re.search(r'Average volume over run is\s+([0-9.eE+-]+)',r.stdout)[1])
    epsilon=1+var*(3.33564e-30)**2/(3*8.8541878128e-12*1.380649e-23*298.15*vol*1e-27)
    gmx_eps=float(re.search(r'Epsilon\s*=\s*([0-9.eE+-]+)',r.stdout)[1])
    assert abs(epsilon-gmx_eps)<5e-4
    data=dict(solvent=meta['solvent'],model='GAFF2/AM1-BCC neutralized to known formal charge',
              last_ps=float(values[-1,0]),analysis_start_ps=a.start_ps,
              density_g_cm3=float(tail[:,columns['Density']].mean()/1000),
              density_temporal_sd_g_cm3=float(tail[:,columns['Density']].std()/1000),
              temperature_K=float(tail[:,columns['Temperature']].mean()),
              epsilon_GROMACS=gmx_eps,epsilon_independent_SI=epsilon,
              blocks_density_g_cm3=[float(y[:,columns['Density']].mean()/1000) for y in np.array_split(tail,4)],
              scope='Single short bulk qualification pilot; no desorption, no experimental efficacy.',
              independent_seeds=1,desorption_model_qualified=False)
    (w/'bulk-analysis.json').write_text(json.dumps(data,indent=2)+'\n');print(data)

if __name__=='__main__':main()
