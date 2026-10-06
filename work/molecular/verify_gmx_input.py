"""Verify actual GROMACS coordinate decoding before any trajectory is used."""
from pathlib import Path
import subprocess,argparse,json
import numpy as np
from ase.io import read
HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory');args=parser.parse_args()
    p=Path(args.directory).resolve()
    r=subprocess.run([str(HERE/'gromacs-env/bin/gmx'),'editconf','-f','initial.gro',
                      '-o','input-parsed.g96'],cwd=p,capture_output=True,text=True,timeout=30)
    (p/'input-decoding.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('GROMACS input parse failed')
    lines=(p/'input-parsed.g96').read_text().splitlines();j=lines.index('POSITION');xyz=[]
    for line in lines[j+1:]:
        if line.strip()=='END':break
        xyz.append([float(x) for x in line.split()[-3:]])
    xyz=np.array(xyz);source=read(p/'packed.xyz').positions*.1
    assert xyz.shape==source.shape
    err=float(np.max(np.abs(xyz-source)))
    assert err<3e-7,err
    d=json.loads((p/'manifest.json').read_text());ns=d['solvent_molecules'];nat=(len(xyz)-d['graphene_atoms']-d['polymer_atoms'])//ns
    sol=read(p/'solvent.xyz').positions*.1
    source_dist=np.linalg.norm(sol[:,None,:]-sol[None,:,:],axis=2)
    block=source[d['graphene_atoms']+d['polymer_atoms']:].reshape(ns,nat,3)
    bond_shape_error=max(float(np.max(np.abs(np.linalg.norm(m[:,None,:]-m[None,:,:],axis=2)-source_dist))) for m in block)
    assert bond_shape_error<5e-7,bond_shape_error
    lines=(p/'index.ndx').read_text().splitlines();groups={};name=None
    for line in lines:
        if line.startswith('['):name=line.strip('[ ]');groups[name]=[]
        elif line.strip():groups[name].extend(int(x) for x in line.split())
    assert len(groups['System'])==len(xyz)
    if d['polymer_n']:
        assert len(groups['GRAPHENE'])==d['graphene_atoms']
        assert len(groups['Mobile'])==len(xyz)-d['graphene_atoms']
    result=dict(coordinates_max_decoding_error_nm=err,
                solvent_rigid_shape_max_difference_nm=bond_shape_error,
                index_group_sizes={k:len(v) for k,v in groups.items()},
                scope='Input serialization and molecular geometry verified, no physical-model qualification')
    (p/'input-verified.json').write_text(json.dumps(result,indent=2)+'\n');print(result)

if __name__=='__main__':main()
