"""Analyse real trajectories: contacts, ester-methyl anchors, H bonding.

No rate, residue-removal success or long-chain extrapolation from a short
single trajectory. Graphene immobility is checked as a model invariant.
"""
from pathlib import Path
import argparse,subprocess,json
import numpy as np
from ase.io import read
from pmma_topology import NB
HERE=Path(__file__).resolve().parent

def frames_g96(path):
    lines=path.read_text().splitlines();time=0.;frames=[];box=None;i=0
    while i<len(lines):
        tag=lines[i].strip();i+=1
        data=[]
        while i<len(lines) and lines[i].strip()!='END':
            data.append(lines[i]);i+=1
        i+=1
        if tag=='TIMESTEP':time=float(data[0].split()[-1])
        elif tag in {'POSITION','POSITIONRED'}:
            xyz=np.array([[float(x) for x in line.split()[-3:]] for line in data])
            frames.append(dict(time_ps=time,xyz_nm=xyz))
        elif tag=='BOX':
            vals=[float(x) for line in data for x in line.split()]
            if len(vals)==3:box=np.diag(vals)
            elif len(vals)==9:
                box=np.array([[vals[0],vals[3],vals[4]],[vals[5],vals[1],vals[6]],[vals[7],vals[8],vals[2]]])
            else:raise RuntimeError('Unexpected box format')
            frames[-1]['box_nm']=box
    return frames

def minimum_image(delta,box):
    f=delta@np.linalg.inv(box)
    return (f-np.rint(f))@box

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory')
    parser.add_argument('--prefix',default='production')
    parser.add_argument('--source',default=None)
    args=parser.parse_args()
    p=Path(args.directory).resolve()
    source=Path(args.source).resolve() if args.source else p
    m=json.loads((source/'manifest.json').read_text())
    ng=m['graphene_atoms'];npma=m['polymer_atoms'];assert ng and npma
    r=subprocess.run([str(HERE/'gromacs-env/bin/gmx'),'trjconv','-s',args.prefix+'.tpr',
                      '-f',args.prefix+'.xtc','-o','analysis-whole.g96','-pbc','whole'],
                     input='0\n',capture_output=True,text=True,cwd=p,timeout=60)
    (p/'trajectory-extraction.txt').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Trajectory extraction failed')
    frames=frames_g96(p/'analysis-whole.g96')
    poly=read(HERE/f'results/pmma/n{m["polymer_n"]}_seed{m["seed"]}/initial.xyz')
    heavy=np.flatnonzero(poly.numbers!=1)
    itp=(HERE/f'results/pmma/n{m["polymer_n"]}_seed{m["seed"]}/pmma.itp').read_text()
    tags=[];flag=False
    for line in itp.splitlines():
        if line.startswith('['):flag=line.strip()=='[ atoms ]';continue
        if flag and line.strip() and not line.startswith(';'):tags.append(line.split()[1][1:])
    assert len(tags)==npma
    methyl=np.array([i for i,t in enumerate(tags) if t=='CE'])
    carbonyl=np.array([i for i,t in enumerate(tags) if t=='O'])
    sigma=np.array([NB[t][1] for t in tags]);eps=np.array([NB[t][2] for t in tags])
    mass=np.array([NB[t][3] for t in tags])
    rows=[];g0=frames[0]['xyz_nm'][:ng].copy();maxmotion=0.
    for f in frames:
        allxyz=f['xyz_nm'];box=f['box_nm'];g=allxyz[:ng];b=allxyz[ng:ng+npma]
        maxmotion=max(maxmotion,float(np.max(np.abs(minimum_image(g-g0,box)))))
        delta=minimum_image(b[:,None,:]-g[None,:,:],box)
        distances=np.linalg.norm(delta,axis=2);nearest=distances.min(axis=1)
        center=np.average(b,axis=0,weights=mass)
        centered=minimum_image(b-center,box)
        rg=float(np.sqrt((mass*np.sum(centered**2,axis=1)).sum()/mass.sum()))
        gyration=(mass[:,None]*centered).T@centered/mass.sum()
        eigenvalues,eigenvectors=np.linalg.eigh(gyration)
        z=float(minimum_image((center-g.mean(axis=0))[None,:],box)[0,2])
        sx=(sigma[:,None]+.347)/2;ex=np.sqrt(eps[:,None]*.275)
        sr6=(sx/distances)**6
        v=4*ex*(sr6**2-sr6)
        vc=4*ex*(sx**12-sx**6) # cutoff1nm, potential-shift default
        epg=float(np.where(distances<1.,v-vc,0.).sum())
        per_atom=np.where(distances<1.,v-vc,0.).sum(axis=1)
        row=dict(time_ps=f['time_ps'],CM_z_from_sheet_nm=z,Rg_nm=rg,
                 Rg_normal_to_sheet_nm=float(np.sqrt(gyration[2,2])),
                 longest_gyration_axis_normal_squared=float(eigenvectors[2,-1]**2),
                 gyration_largest_minus_middle_nm2=float(eigenvalues[-1]-eigenvalues[-2]),
                 nearest_heavy_min_nm=float(nearest[heavy].min()),
                 heavy_contacts_within045nm=int(np.sum(nearest[heavy]<.45)),
                 ester_methyl_contacts_within045nm=int(np.sum(nearest[methyl]<.45)),
                 smooth_heavy_contacts=float(np.sum(1/(1+(nearest[heavy]/.45)**6))),
                 polymer_graphene_LJ_energy_kJ_mol=epg)
        for group,accepted in [('ester_methyl',{'CE','H1'}),
                               ('ester_core',{'C','O','OS'}),
                               ('backbone_and_side_methyl',{'CX','CA','HC'})]:
            row[group+'_G_LJ_kJ_mol']=float(sum(per_atom[i] for i,t in enumerate(tags) if t in accepted))
        if m['solvent']=='HFIP':
            solv=allxyz[ng+npma:].reshape(m['solvent_molecules'],12,3)
            h=solv[:,10];donor=solv[:,0];acceptor=b[carbonyl]
            ha=minimum_image(acceptor[None,:,:]-h[:,None,:],box)
            hd=minimum_image(donor-h,box)[:,None,:]
            dr=np.linalg.norm(ha,axis=2)
            cosine=np.sum(ha*hd,axis=2)/(np.linalg.norm(hd,axis=2)*dr)
            row['HFIP_carbonyl_hbonds']=int(np.sum((dr<.25)&(cosine<-.8660254038)))
        rows.append(row)
    # XTC precision is10000/nm. Periodic wrapping before quantization can
    # change the apparent coordinate by ~1e-4nm even for a frozen atom.
    # This checks immobility only to the saved trajectory resolution.
    assert maxmotion<1.1e-4,maxmotion
    tail=[r for r in rows if r['time_ps']>=max(0,rows[-1]['time_ps']/2)]
    numeric=[k for k,v in rows[0].items() if isinstance(v,(int,float)) and k!='time_ps']
    stats={k:dict(mean=float(np.mean([r[k] for r in tail])),std=float(np.std([r[k] for r in tail]))) for k in numeric}
    result=dict(scope='Short exploratory oligomer/rigid-sheet trajectory. Not950K, not cleaning, not oxide safety, not a converged desorption free energy.',
                solvent=m['solvent'],polymer_n=m['polymer_n'],frames=len(rows),
                biased_trajectory=(args.prefix!='production'),
                last_ps=rows[-1]['time_ps'],graphene_max_apparent_motion_nm=maxmotion,
                saved_coordinate_resolution_nm=1e-4,
                rigidity_check_is_not_material_damage_prediction=True,
                statistics_second_half=stats,rows=rows,independent_replicas=1)
    (p/'interface-analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(m['solvent'],rows[-1]['time_ps'],stats,'Gmotion',maxmotion)

if __name__=='__main__':main()
