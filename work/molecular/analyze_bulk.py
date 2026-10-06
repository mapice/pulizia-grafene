"""Read GROMACS column legends, verify dielectric fluctuation formula in SI.

No inference of convergence from trajectory length alone. All short runs are
pilots; a solvent-density agreement cannot qualify interface selectivity.
"""
import subprocess,json,re,argparse
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
P=HERE/'results/gromacs/bulk_scale0.5_seed20261006'
GMX=str(HERE/'gromacs-env/bin/gmx')

def xvg(path):
    columns={};rows=[]
    for line in path.read_text().splitlines():
        m=re.search(r'@\s*s(\d+)\s+legend\s+"([^"]+)"',line)
        if m:columns[m[2]]=int(m[1])+1
        if line.strip() and not line.startswith(('@','#')):rows.append([float(x) for x in line.split()])
    return np.array(rows),columns

def main():
    global P
    parser=argparse.ArgumentParser()
    parser.add_argument('--scale14',type=float,default=.5)
    parser.add_argument('--seed',type=int,default=20261006)
    parser.add_argument('--start-ps',type=float,default=200)
    args=parser.parse_args()
    P=HERE/'results/gromacs'/f'bulk_scale{args.scale14:g}_seed{args.seed}'
    e=subprocess.run([GMX,'energy','-f','npt.edr','-o','thermo.xvg'],input='Temperature\nDensity\nVolume\n0\n',text=True,capture_output=True,cwd=P,timeout=30)
    (P/'thermo-analysis.txt').write_text(e.stdout+e.stderr)
    if e.returncode:raise RuntimeError('Energy analysis failed')
    dip=subprocess.run([GMX,'dipoles','-s','npt.tpr','-f','npt.xtc','-b',str(args.start_ps),'-o','dipole.xvg','-eps','epsilon.xvg','-temp','298.15','-nopairs'],input='0\n',capture_output=True,text=True,cwd=P,timeout=30)
    (P/'dipole-analysis.txt').write_text(dip.stdout+dip.stderr)
    if dip.returncode:raise RuntimeError('Dipole analysis failed')
    a,c=xvg(P/'thermo.xvg');y=a[a[:,0]>=args.start_ps]
    assert all(k in c for k in ['Temperature','Density','Volume']),c
    rho=y[:,c['Density']]/1000;t=y[:,c['Temperature']];v=y[:,c['Volume']]
    d,dc=xvg(P/'dipole.xvg');m=d[:,1:4]
    var=float((m*m).sum(axis=1).mean()-(m.mean(axis=0)**2).sum())
    # Use exactly GROMACS average volume over dipole frames for comparison.
    vv=float(re.search(r'Average volume over run is\s+([0-9.eE+-]+)',dip.stdout)[1])
    eps=1+var*(3.33564e-30)**2/(3*8.8541878128e-12*1.380649e-23*298.15*vv*1e-27)
    gmx_eps=float(re.search(r'Epsilon\s*=\s*([0-9.eE+-]+)',dip.stdout)[1])
    assert abs(eps-gmx_eps)<5e-5,(eps,gmx_eps)
    blocks=[]
    for part in np.array_split(y,4):blocks.append(dict(time_start_ps=float(part[0,0]),time_end_ps=float(part[-1,0]),density_g_cm3=float(part[:,c['Density']].mean()/1000),temperature_K=float(part[:,c['Temperature']].mean())))
    r=dict(scope='Bulk HFIP pilot, not interface cleaning. One seed, 125 molecules; dielectric not a converged experimental prediction.',scale14=args.scale14,seed=args.seed,last_ps=float(a[-1,0]),analysis_start_ps=args.start_ps,density_g_cm3=float(rho.mean()),density_sd_g_cm3=float(rho.std()),temperature_K=float(t.mean()),volume_nm3=float(v.mean()),dielectric_GROMACS=gmx_eps,dielectric_independent_SI=eps,dipole_variance_Debye2=var,dipole_frames=len(d),density_reference_experimental_g_cm3=1.607,dielectric_reference_experimental=16.7,reference_doi='10.1002/psc.3543',blocks=blocks,forcefield_qualified_for_desorption=False,qualification_reason='Density alone insufficient; response, convergence, conventions and interfacial interactions require validation before quantitative solvent ranking.')
    (P/'summary.json').write_text(json.dumps(r,indent=2)+'\n');print(r)

if __name__=='__main__':main()
