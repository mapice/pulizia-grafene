"""Executable second sampling coordinate from an observed hidden mode.

Mass-weighted normal gyration radius, rather than a discontinuous atom
contact or eigenvector. Restricted here to the already checked n4 chain.
Validate values against independent trajectory analysis and bias forces
against numerical derivatives before a biased dynamics pilot.
"""
from pathlib import Path
import argparse,json,subprocess,os
import numpy as np
from analyze_interface import frames_g96

HERE=Path(__file__).resolve().parent
BASE=HERE/'results/gaff-systems/HFIP_n4_seed20261006'


def masses():
    data=[];flag=False
    for line in (HERE/'results/pmma/n4_seed20261006/pmma.itp').read_text().splitlines():
        if line.startswith('['):flag=line.strip()=='[ atoms ]';continue
        if flag and line.strip() and not line.startswith(';'):data.append(float(line.split()[7]))
    assert len(data)==65
    return np.array(data)


def specification(file='CV',stride=1,bias=None):
    mass=masses();ng=json.loads((BASE/'manifest.json').read_text())['graphene_atoms']
    indices=list(range(ng+1,ng+len(mass)+1))
    weights=','.join(f'{m:.12g}' for m in mass)
    lines=[f'WHOLEMOLECULES ENTITY0={ng+1}-{ng+len(mass)}',
           f'g: CENTER ATOMS=1-{ng} NOPBC',
           f'p: CENTER ATOMS={ng+1}-{ng+len(mass)} WEIGHTS={weights} NOPBC',
           'height: DISTANCE ATOMS=g,p COMPONENTS']
    for i in indices:lines.append(f'd{i}: DISTANCE ATOMS=p,{i} COMPONENTS NOPBC')
    args=','.join(f'd{i}.z' for i in indices)
    coeff=','.join(f'{m/mass.sum():.16g}' for m in mass)
    lines += [f'varz: COMBINE ARG={args} COEFFICIENTS={coeff} POWERS='+','.join(['2']*len(indices))+' PERIODIC=NO',
              'rgz: CUSTOM ARG=varz FUNC=sqrt(x) PERIODIC=NO']
    if bias:
        z,rg,kz,krg=bias
        lines.append(f'bias: RESTRAINT ARG=height.z,rgz AT={z},{rg} KAPPA={kz},{krg}')
    lines.append(f'PRINT ARG=height.z,rgz'+(',bias.bias' if bias else '')+f' STRIDE={stride} FILE={file} FMT=%16.10f')
    return '\n'.join(lines)+'\n'


def main():
    p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');args=p.parse_args()
    root=HERE/'results/orientation-control';root.mkdir(exist_ok=True)
    (root/'monitor.dat').write_text(specification())
    if not args.validate:return
    rows=[]
    env=os.environ.copy();env['OMP_NUM_THREADS']='1'
    for direction in ['forward','reverse']:
        folder=HERE/'results/umbrella-ladders/HFIP_n4_k750_500ps'/direction/'z1.050'
        assert (folder/'complete.json').exists()
        cmd=[str(HERE/'gromacs-env/bin/plumed'),'driver','--ixtc',str(folder/'run.xtc'),
             '--plumed','monitor.dat','--timestep','.001','--trajectory-stride','1000']
        r=subprocess.run(cmd,cwd=root,env=env,capture_output=True,text=True,timeout=60)
        (root/f'{direction}-driver.txt').write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError('PLUMED monitor failed')
        x=np.loadtxt(root/'CV',comments='#')
        target=root/f'{direction}-CV.dat';(root/'CV').replace(target)
        independent=json.loads((folder/'interface-analysis.json').read_text())['rows']
        assert len(x)==len(independent)
        refs=np.array([[r['time_ps'],r['CM_z_from_sheet_nm'],r['Rg_normal_to_sheet_nm']] for r in independent])
        errors=np.max(np.abs(x-refs),axis=0)
        assert errors[0]<1e-6 and max(errors[1:])<3e-6,errors
        rows.append(dict(direction=direction,max_errors_time_z_rgz=errors.tolist(),frames=len(x)))
    # One whole frame, retaining only surface/polymer; solvent atoms do not
    # enter these two coordinates or their explicit restraint derivatives.
    frame=frames_g96(HERE/'results/umbrella-ladders/HFIP_n4_k750_500ps/forward/z1.050/analysis-whole.g96')[0]
    xyz=frame['xyz_nm'][:513]
    text=f'{len(xyz)}\nOne frame, coordinates nm; not a dynamics run\n'
    text+='\n'.join('C '+ ' '.join(f'{v:.10f}' for v in atom) for atom in xyz)+'\n'
    (root/'one-frame.xyz').write_text(text)
    (root/'bias-force-check.dat').write_text(specification(file='bias-CV',bias=(1.05,.15,750,2000)))
    cmd=[str(HERE/'gromacs-env/bin/plumed'),'driver','--ixyz','one-frame.xyz',
         '--plumed','bias-force-check.dat','--box',','.join(map(str,np.diag(frame['box_nm']))),
         '--debug-forces','bias-force-derivatives.dat']
    r=subprocess.run(cmd,cwd=root,env=env,capture_output=True,text=True,timeout=120)
    (root/'force-driver.txt').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('PLUMED derivative check failed')
    derivatives=np.loadtxt(root/'bias-force-derivatives.dat',comments='#')
    finite_error=float(np.abs(derivatives[:,1]-derivatives[:,2]).max())
    assert finite_error<.001,finite_error
    mass=masses();wm=mass/mass.sum();poly=xyz[448:]
    center=np.average(poly,axis=0,weights=mass)
    height=center[2]-xyz[:448,2].mean()
    delta=poly[:,2]-center[2];rg=float(np.sqrt(np.sum(wm*delta**2)))
    independent=np.zeros((513,3))
    independent[:448,2]=750*(height-1.05)/448
    independent[448:,2]=-750*(height-1.05)*wm-2000*(rg-.15)*wm*delta/rg
    analytical_error=float(np.abs(independent.ravel()-derivatives[:1539,1]).max())
    assert analytical_error<6e-7,analytical_error
    assert np.abs(independent.sum(axis=0)).max()<1e-10
    result=dict(coordinate_value_checks=rows,
                derivative_file='bias-force-derivatives.dat',
                finite_difference_max_error_kJ_mol_nm=finite_error,
                independent_analytical_force_error_kJ_mol_nm=analytical_error,
                analytical_forces_net_sum_zero=True,
                coordinate_and_first_derivative_validation_passed=True,
                scope='Coordinate and force implementation checks, not desorption equilibrium or cleaning',
                qualified_for_cleaning=False)
    (root/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__':main()
