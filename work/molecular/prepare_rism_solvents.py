"""Independent statistical-liquid route: DRISM with measured dielectric input.

MDL units verified in Amber26 manual§8.7.1: charge*18.2223, epsilon kcal,
rmin/2 Angstrom, coordinates Angstrom. Exactly equivalent charge/LJ/mass
types share a site; physical coordinates are reordered by site group. Rigid
conformer is an explicit approximation; dielectric is imposed experimental
input. No organic susceptibility has yet met the convergence test.
"""
from pathlib import Path
import json,time,argparse,subprocess,os
import numpy as np
import parmed
from parmed.constants import AMBER_ELECTROSTATIC

HERE=Path(__file__).resolve().parent
AMBER=HERE/'amber-env'
PROP={'HFIP':(1.607,16.7,'10.1002/psc.3543'),
      'DMF':(.944,37.2,'10.1021/je9010773')}

def field(name,values,kind):
    if kind=='int':
        fmt='%FORMAT(10I8)';lines=[''.join(f'{int(x):8d}' for x in values[j:j+10]) for j in range(0,len(values),10)]
    elif kind=='char':
        fmt='%FORMAT(20a4)';lines=[''.join(f'{str(x):<4s}'[:4] for x in values[j:j+20]) for j in range(0,len(values),20)]
    else:
        fmt='%FORMAT(5e16.8)';lines=[''.join(f'{float(x):16.8E}' for x in values[j:j+5]) for j in range(0,len(values),5)]
    return '%FLAG '+name+'\n'+fmt+'\n'+'\n'.join(lines)+'\n'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--solvent',choices=PROP,default='HFIP');a=parser.parse_args()
    name=a.solvent;rho,dieps,ref=PROP[name];src=HERE/'results/gaff2'/name
    s=parmed.load_file(str(src/'molecule.prmtop'),xyz=str(src/'molecule.inpcrd'))
    n=len(s.atoms);mass=sum(x.mass for x in s.atoms);charge=sum(x.charge for x in s.atoms)
    equivalent={}
    for i,x in enumerate(s.atoms):
        key=(round(x.charge,8),round(x.epsilon,8),round(x.rmin,8),x.mass)
        equivalent.setdefault(key,[]).append(i)
    groups=list(equivalent.values());representatives=[s.atoms[g[0]] for g in groups]
    ns=len(groups)
    # Standard MDL order groups physical coordinates by equivalent site type.
    coords=np.concatenate([s.coordinates[g] for g in groups],axis=0)
    assert abs(charge)<1e-7
    p=HERE/'results/rism'/name;p.mkdir(parents=True,exist_ok=True)
    mdl='%VERSION  VERSION_STAMP = V0001.000 DATE = 10/06/26 00:00:00\n'
    mdl+=field('TITLE',[name],'char')
    mdl+=field('POINTERS',[n,ns],'int')+field('ATMTYP',range(1,ns+1),'int')
    mdl+=field('ATMNAME',[x.name for x in representatives],'char')
    mdl+=field('MASS',[x.mass for x in representatives],'real')
    mdl+=field('CHG',[x.charge*AMBER_ELECTROSTATIC for x in representatives],'real')
    mdl+=field('LJEPSILON',[x.epsilon for x in representatives],'real')
    mdl+=field('LJSIGMA',[x.rmin for x in representatives],'real')
    mdl+=field('MULTI',[len(g) for g in groups],'int')+field('COORD',coords.ravel(),'real')
    (p/'solvent.mdl').write_text(mdl)
    molarity=1000*rho/mass
    inp=f'''&PARAMETERS
  THEORY='DRISM', CLOSURE='KH',
  NR=8192, DR=0.025,
  OUTLIST='xgct', ROUT=0, KOUT=0,
  MDIIS_NVEC=20, MDIIS_DEL=0.3, TOLERANCE=1.e-10,
  KSAVE=-1, PROGRESS=100, MAXSTEP=10000,
  SMEAR=1, ADBCOR=0.5,
  TEMPERATURE=298.15, DIEPS={dieps}, NSP=1,
  entropicDecomp=0
/
&SPECIES
  DENSITY={molarity:.12f}, UNITS='M',
  MODEL='solvent.mdl'
/
'''
    (p/'liquid.inp').write_text(inp)
    manifest=dict(solvent=name,method='DRISM/KH',sites=ns,physical_atoms=n,
                  equivalent_site_multiplicities=[len(g) for g in groups],
                  rigid_conformer=True,geometry_source='neutral GAFF2 input geometry',
                  density_g_cm3=rho,molarity_mol_L=molarity,
                  dielectric_imposed_from_experiment=dieps,reference_doi=ref,
                  dielectric_is_not_predicted=True,
                  mdl_LJ_radius='rmin/2 Angstrom, not sigma',
                  mdl_charge_factor=AMBER_ELECTROSTATIC,
                  scope='Prepare bulk susceptibility for independent solvation-theory check; no cleaning result')
    (p/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(p,manifest,flush=True)

if __name__=='__main__':main()
