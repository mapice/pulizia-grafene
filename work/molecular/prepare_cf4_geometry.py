"""GFN2 geometry subprocess; do not load conflicting torch/OpenMP here."""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
from pathlib import Path
import sys,json
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import FIRE
from tblite.ase import TBLite
root=Path(sys.argv[1]);assert root.is_dir()
vertices=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],dtype=float)/np.sqrt(3)
a=Atoms('CF4',positions=np.vstack([np.zeros(3),vertices*1.33]))
a.calc=TBLite(method='GFN2-xTB',verbosity=0)
converged=bool(FIRE(a,logfile=str(root/'geometry-relaxation.txt'),maxstep=.05).run(fmax=.005,steps=120))
assert converged
fmax=float(np.linalg.norm(a.get_forces(),axis=1).max())
a.positions-=a.positions[0];a.calc=None;write(root/'monomer.xyz',a)
(root/'geometry.json').write_text(json.dumps(dict(GFN2_geometry_converged=converged,
    fmax_eV_A=fmax,starting_bond_length_A=1.33,
    relaxed_CF_bond_lengths_A=np.linalg.norm(a.positions[1:],axis=1).tolist()),indent=2)+'\n')
