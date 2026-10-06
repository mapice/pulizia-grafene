"""Explicit transcription: Marchelli et al. 2022 SI Tables I-V.

Lorentz-Berthelot LJ, kcal/Angstrom converted to kJ/nm. LAMMPS harmonic
coefficients lack the OpenMM 1/2 prefactor, so k is doubled. Intramolecular
1-4 scaling was not specified in the article: it is an explicit sensitivity
parameter, never silently declared verified. No interfacial qualification.
"""
import numpy as np
from rdkit import Chem
import openmm as mm
from openmm import unit as u

SMILES='OC(C(F)(F)F)C(F)(F)F'
MOL=Chem.AddHs(Chem.MolFromSmiles(SMILES))
TYPES=[]
for a in MOL.GetAtoms():
    z=a.GetSymbol()
    TYPES.append('OH' if z=='O' else 'HO' if z=='H' and a.GetNeighbors()[0].GetSymbol()=='O' else 'HC' if z=='H' else 'F' if z=='F' else 'CF' if sum(n.GetSymbol()=='F' for n in a.GetNeighbors())==3 else 'CT')
# charge [e], sigma [nm], epsilon [kJ/mol]
NONBONDED={'CF':(.6,.33611,.097098*4.184),'CT':(-.07,.33611,.097098*4.184),'F':(-.2,.31578,.071041*4.184),'OH':(-.595,.29548,.203255*4.184),'HO':(.495,0.,0.),'HC':(.17,.23734,.028321*4.184)}
BONDS={tuple(sorted(k)):(r*.1,2*kr*4.184*100) for k,kr,r in [(('CT','F'),734,1.36),(('CT','CT'),268,1.53),(('CT','OH'),320,1.36),(('CT','HC'),2000,1.09),(('OH','HO'),2000,1.0)]}
def canonical(t):return min(tuple(t),tuple(t)[::-1])
ANGLES={canonical(k):(np.deg2rad(theta),2*kr*4.184) for k,kr,theta in [(('F','CT','F'),55.048,107.6),(('F','CT','CT'),55.048,111),(('CT','CT','CT'),55.048,110),(('CT','CT','OH'),55.048,111),(('CT','CT','HC'),55.048,109.5),(('OH','CT','HC'),47.548,109.5),(('CT','OH','HO'),47.548,109.5)]}
MASSES={'H':1.008,'C':12.011,'O':15.9994,'F':18.9984}
def bonded_type(t):return 'CT' if t=='CF' else t

def system(nmolecules,box,scale14=.5,cutoff=1.0):
    s=mm.System()
    s.setDefaultPeriodicBoxVectors(mm.Vec3(box,0,0),mm.Vec3(0,box,0),mm.Vec3(0,0,box))
    nb=mm.NonbondedForce();nb.setNonbondedMethod(nb.PME);nb.setCutoffDistance(cutoff*u.nanometer);nb.setEwaldErrorTolerance(1e-4)
    nb.setUseDispersionCorrection(True)
    bf=mm.HarmonicBondForce();af=mm.HarmonicAngleForce()
    allbonds=[]
    for m in range(nmolecules):
        offset=m*MOL.GetNumAtoms()
        for a,t in zip(MOL.GetAtoms(),TYPES):
            s.addParticle(MASSES[a.GetSymbol()]);nb.addParticle(*NONBONDED[t])
        for b in MOL.GetBonds():
            i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            key=tuple(sorted((bonded_type(TYPES[i]),bonded_type(TYPES[j]))))
            r,k=BONDS[key]
            if MOL.GetAtomWithIdx(i).GetSymbol()=='H' or MOL.GetAtomWithIdx(j).GetSymbol()=='H':
                s.addConstraint(offset+i,offset+j,r)
            else:bf.addBond(offset+i,offset+j,r,k)
            allbonds.append((offset+i,offset+j))
        for a in MOL.GetAtoms():
            ns=[n.GetIdx() for n in a.GetNeighbors()]
            for ii,i in enumerate(ns):
                for k in ns[ii+1:]:
                    j=a.GetIdx();key=canonical([bonded_type(TYPES[x]) for x in [i,j,k]])
                    theta,force=ANGLES[key]
                    af.addAngle(offset+i,offset+j,offset+k,theta,force)
    nb.createExceptionsFromBonds(allbonds,scale14,scale14)
    s.addForce(nb);s.addForce(bf);s.addForce(af)
    # HFIP dihedral coefficients are zero in the published model.
    if abs(sum(NONBONDED[t][0] for t in TYPES))>1e-10:raise ValueError('Nonneutral molecule')
    return s,np.array([NONBONDED[t][0] for t in TYPES]),sum(MASSES[a.GetSymbol()] for a in MOL.GetAtoms())
