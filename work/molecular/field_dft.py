"""Independent fixed-nuclei electric response of the actual HFIP proxy.

PBE0/def2-TZVP at the same geometry used in learned-model checks. Uniform
field enters the electronic one-body Hamiltonian; nuclear field energy is
added explicitly. D3 is field-independent at fixed nuclei and cancels in
the curvature. This is gas-phase response, not bulk permittivity.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='2'
os.nice(8)
from pathlib import Path
import time,json,resource,sys
import numpy as np
from ase.io import read
from pyscf import gto,dft,lib
lib.num_threads(2)
HERE=Path(__file__).resolve().parent;OUT=HERE/'results/field-dft';OUT.mkdir(parents=True,exist_ok=True)
HARTREE_EV=27.211386245988
FIELD_AU_V_A=51.4220674763
a=read(HERE/'results/clusters/HFIP.refined.xyz')

def main():
    rows=[];dm=None;start=time.monotonic()
    for field in [0.,-.002,.002,-.01,.01]:
        mol=gto.M(atom=list(zip(a.get_chemical_symbols(),a.positions)),basis='def2-tzvp',unit='Angstrom',spin=0,charge=0,max_memory=1500,verbose=4,output=str(OUT/f'field{field:+.3f}.scf.log'))
        mf=dft.RKS(mol).density_fit();mf.xc='pbe0';mf.grids.level=3;mf.conv_tol=1e-11;mf.max_cycle=100;mf.chkfile=None
        base_h=mf.get_hcore();rz=mol.intor('int1e_r',comp=3)[2];au_field=field/FIELD_AU_V_A
        mf.get_hcore=lambda *args,**kwargs:base_h+au_field*rz
        e=float(mf.kernel(dm0=dm));dm=mf.make_rdm1()
        if not mf.converged:raise RuntimeError('Field SCF unconverged')
        enuclear_field=-au_field*float(np.dot(mol.atom_charges(),mol.atom_coords()[:,2]))
        mu=np.array(mf.dip_moment(unit='Debye',verbose=0))
        row=dict(field_V_A=field,total_energy_eV=(e+enuclear_field)*HARTREE_EV,dipole_Debye=mu.tolist(),SCF_converged=bool(mf.converged),SCF_energy_Ha=e,nuclear_field_Ha=enuclear_field)
        rows.append(row);print(row,flush=True)
    result=[]
    for j,delta in [(1,.002),(3,.01)]:
        curvature=-(rows[j]['total_energy_eV']+rows[j+1]['total_energy_eV']-2*rows[0]['total_energy_eV'])/delta**2
        linear=(rows[j+1]['total_energy_eV']-rows[j]['total_energy_eV'])/(2*delta)
        result.append(dict(field_magnitude_V_A=delta,polarizability_from_energy_eA2_V=curvature,energy_linear_derivative_eA=linear,dipole_zero_eA=rows[0]['dipole_Debye'][2]/4.8032045098))
        assert curvature>0,'Stable electronic response must be positive'
    data=dict(method='PBE0/def2-TZVP',scope='Fixed-nuclei gas HFIP field response; not liquid dielectric',geometry='HFIP.refined.xyz',rows=rows,response=result,elapsed_seconds=time.monotonic()-start,peak_RSS_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024))
    (OUT/'response.json').write_text(json.dumps(data,indent=2)+'\n');print(result,flush=True)

if __name__=='__main__':main()
