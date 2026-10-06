"""Frozen-geometry HFIP OH torsion: quantum and published classical models.

Only the OH hydrogen moves. Bond length and C-O-H angle remain fixed, so
classical bonded terms are constant and cancel in relative energy. This is
a diagnostic gas rotor; no liquid free energy or cleaning prediction.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[v]='2'
os.nice(8)
from pathlib import Path
import time,json,resource,sys,hashlib
import numpy as np
from ase.io import read,write
from rdkit import Chem
from pyscf import gto,dft,lib
from dftd3.interface import DispersionModel,RationalDampingParam
from tblite.ase import TBLite
import openmm as mm
from openmm import unit as u
from published_hfip import MOL,TYPES,NONBONDED,system

lib.num_threads(2)
HERE=Path(__file__).resolve().parent
OUT=HERE/'results/hfip-rotor';OUT.mkdir(parents=True,exist_ok=True)
HA_KJ=2625.499638;EV_KJ=96.4853321233;BOHR_A=.529177210903
COULOMB=138.93545764438198
ho=TYPES.index('HO');oh=TYPES.index('OH');hc=TYPES.index('HC');ct=TYPES.index('CT')
original=read(HERE/'results/clusters/HFIP.refined.xyz')
dist=Chem.GetDistanceMatrix(MOL)
params=np.array([NONBONDED[t] for t in TYPES])

def geometry(angle):
    a=original.copy()
    # Reverse the conventional H_O-O-C-H_C quartet so the moved OH-H is last.
    a.set_dihedral(hc,ct,oh,ho,float(angle))
    actual=a.get_dihedral(ho,oh,ct,hc)
    assert abs((actual-angle+180)%360-180)<1e-7,(actual,angle)
    assert abs(a.get_distance(ho,oh)-original.get_distance(ho,oh))<1e-12
    assert abs(a.get_angle(ct,oh,ho)-original.get_angle(ct,oh,ho))<1e-10
    return a

def classical(a,scale):
    value=0.
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            if dist[i,j]<=2:continue
            factor=scale if dist[i,j]==3 else 1.
            r=np.linalg.norm(a.positions[i]-a.positions[j])*.1
            q1,s1,e1=params[i];q2,s2,e2=params[j]
            sr6=((s1+s2)/2/r)**6
            value+=factor*(4*np.sqrt(e1*e2)*(sr6*sr6-sr6)+COULOMB*q1*q2/r)
    return float(value)

def dft_energy(a,angle,dm=None):
    stem=f'angle{int(angle):+04d}'
    digest=hashlib.sha256(a.positions.tobytes()).hexdigest()
    path=OUT/(stem+'.json')
    if path.exists():
        r=json.loads(path.read_text())
        if r['geometry_sha256']==digest and r['converged']:
            return r,None
    mol=gto.M(atom=list(zip(a.get_chemical_symbols(),a.positions)),basis='def2-tzvp',
              unit='Angstrom',charge=0,spin=0,verbose=4,max_memory=1500,
              output=str(OUT/(stem+'.scf.log')))
    mf=dft.RKS(mol).density_fit()
    mf.xc='pbe0';mf.grids.level=3;mf.conv_tol=1e-10;mf.max_cycle=100;mf.chkfile=None
    start=time.monotonic()
    energy=float(mf.kernel(dm0=dm))
    if not mf.converged:raise RuntimeError('Unconverged torsion '+stem)
    correction=float(DispersionModel(a.numbers,a.positions/BOHR_A).get_dispersion(
        RationalDampingParam(method='pbe0'),grad=False)['energy'])
    mu=np.array(mf.dip_moment(unit='Debye',verbose=0))
    r=dict(angle_degree=float(angle),method='PBE0-D3(BJ)/def2-TZVP',
           electronic_Ha=energy,dispersion_Ha=correction,
           total_kJ_mol=(energy+correction)*HA_KJ,dipole_Debye=mu.tolist(),
           dipole_norm_Debye=float(np.linalg.norm(mu)),converged=True,
           geometry_sha256=digest,elapsed_seconds=time.monotonic()-start)
    path.write_text(json.dumps(r,indent=2)+'\n')
    return r,mf.make_rdm1()

def main():
    rows=[];start=time.monotonic();dm=None
    contexts={}
    for scale in [.5,1.]:
        s,_,_=system(1,10.,scale14=scale)
        nb=s.getForce(0);nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        contexts[scale]=mm.Context(s,mm.VerletIntegrator(.001),
                                  mm.Platform.getPlatformByName('Reference'))
    for angle in range(-180,180,30):
        a=geometry(angle);write(OUT/f'angle{angle:+04d}.xyz',a)
        row,dm=dft_energy(a,angle,dm)
        a.calc=TBLite(method='GFN2-xTB',verbosity=0,accuracy=1.0)
        row['GFN2_kJ_mol']=float(a.get_potential_energy())*EV_KJ
        for scale,context in contexts.items():
            row[f'classical{scale:g}_nonbonded_kJ_mol']=classical(a,scale)
            context.setPositions(a.positions*.1*u.nanometer)
            row[f'OpenMM{scale:g}_total_kJ_mol']=float(context.getState(
                getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole))
        rows.append(row)
        print(angle,row['dipole_norm_Debye'],row['elapsed_seconds'],flush=True)
        (OUT/'partial.json').write_text(json.dumps(rows,indent=2)+'\n')
    fields=['total_kJ_mol','GFN2_kJ_mol',
            'classical0.5_nonbonded_kJ_mol','classical1_nonbonded_kJ_mol']
    for field in fields:
        minimum=min(r[field] for r in rows)
        for row in rows:row[field.replace('_kJ_mol','_relative_kJ_mol')]=row[field]-minimum
    errors={}
    for scale in [.5,1.]:
        dclass=np.array([r[f'classical{scale:g}_nonbonded_kJ_mol'] for r in rows])
        dmm=np.array([r[f'OpenMM{scale:g}_total_kJ_mol'] for r in rows])
        err=float(np.max(np.abs((dclass-dclass[0])-(dmm-dmm[0]))))
        assert err<1e-8,err
        errors[str(scale)]=err
    result=dict(scope='Gas HFIP frozen OH rotor; not liquid free energy or cleaning.',
                quartet_zero_based=[ho,oh,ct,hc],initial_angle_degree=original.get_dihedral(ho,oh,ct,hc),
                fixed_nuclei_except='OH-H only',rows=rows,
                direct_vs_OpenMM_relative_energy_max_error_kJ_mol=errors,
                elapsed_seconds=time.monotonic()-start,
                peak_RSS_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/
                    (1024**2 if sys.platform=='darwin' else 1024))
    (OUT/'scan.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Finished',result['elapsed_seconds'],result['peak_RSS_MB'],errors,flush=True)

if __name__=='__main__':main()
