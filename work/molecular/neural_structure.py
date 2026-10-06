"""Shared immutable input reconstruction, with no thread/process side effects."""
from pathlib import Path
import json
import numpy as np
from ase import Atoms
from ase.io import read
HERE=Path(__file__).resolve().parent
SHA='fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'


def structure(name):
    if name=='monomer':
        a=read(HERE/'results/clusters/HFIP.refined.xyz');a.positions-=a.positions.mean(0);return a
    assert name in ['liquid8','liquid64']
    p=HERE/'results/gaff-systems/HFIP_n0_seed20261008_N64'
    status=json.loads((p/'run-status.json').read_text())
    assert status[-1]['stage']=='production' and status[-1]['exitcode']==0
    mono=read(HERE/'results/clusters/HFIP.refined.xyz');gro=read(p/'production.gro')
    count=8 if name=='liquid8' else 64
    a=Atoms(numbers=np.tile(mono.numbers,count),positions=gro.positions[:count*12],cell=gro.cell,pbc=True)
    assert len(a)==12*count and min(a.cell.lengths())>20
    return a
