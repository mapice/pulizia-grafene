"""Independent double-precision cell-derivative audit of the polar model.

Keep forces/energy fixed, compare stress with central finite strains at
two step sizes. Exact scalar energy comparisons precede NPT acceptance.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.set_default_dtype(torch.float64)
from ase.io import read
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',required=True);args=p.parse_args()
    root=HERE/'results/virial-geometry-audit';root.mkdir(exist_ok=True)
    target=root/f'{args.tag}.json';assert not target.exists()
    source=HERE/'results/neural-bulk/HFIP_N16_seed20261016/dynamics-pilot/relaxed-initial.xyz'
    a=read(source);a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).double()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')
    a.calc=calc;start=time.monotonic()
    e=float(a.get_potential_energy());forces=a.get_forces().copy();stress=a.get_stress().copy();volume=a.get_volume()
    checks=[]
    for kind in ['xx','isotropic']:
        for delta in [1e-4,5e-5]:
            values=[]
            for sign in [-1,1]:
                b=a.copy();cell=a.cell.array.copy()
                if kind=='xx':cell[0,:]*=1+sign*delta
                else:cell*=1+sign*delta
                b.set_cell(cell,scale_atoms=True);b.calc=calc
                values.append(float(b.get_potential_energy()))
            numeric=(values[1]-values[0])/(2*delta*volume*(1 if kind=='xx' else 3))
            analytic=float(stress[0] if kind=='xx' else stress[:3].mean())
            row=dict(strain=kind,step=delta,analytic_eV_A3=analytic,numerical_eV_A3=numeric,
                     difference_eV_A3=numeric-analytic,energies_eV=values)
            checks.append(row);print(row,flush=True)
    result=dict(model_sha256=SHA,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                extension_sha256=hashlib.sha256((HERE/'.venv/lib/python3.12/site-packages/mace/modules/extensions.py').read_bytes()).hexdigest(),
                dtype='float64',atoms=len(a),energy_eV=e,forces_eV_A=forces.tolist(),stress_eV_A3=stress.tolist(),
                cell_A=a.cell.array.tolist(),volume_A3=float(volume),checks=checks,
                elapsed_seconds=time.monotonic()-start,scope='Implementation audit, not liquid-property or cleaning qualification')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('Audit completed',args.tag,'max residual',max(abs(c['difference_eV_A3']) for c in checks),flush=True)


if __name__=='__main__':main()
