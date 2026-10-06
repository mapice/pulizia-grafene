"""One additional same-model gas pair, retaining previous comparisons."""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.set_default_dtype(torch.float64)
from ase.io import read
from mace.calculators import MACECalculator
from neural_structure import SHA
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--pair',required=True);args=p.parse_args()
    root=HERE/'results/polar-size-audit';target=root/(args.pair+'-M.json');assert not target.exists()
    record=json.loads((HERE/'results/clusters/binding-refined.json').read_text())['pairs'][args.pair]
    a=read(HERE/'results/clusters'/record['best_geometry'])
    n=next(t['n_acceptor'] for t in record['trials'] if t['name']+'.refined.xyz'==record['best_geometry'])
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).double()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')
    def energy(b):
        b=b.copy();b.positions-=b.positions.mean(0);b.info.update(charge=0,spin=1,external_field=[0.,0.,0.]);b.calc=calc
        return float(b.get_potential_energy())
    ab,fa,fb=energy(a),energy(a[:n]),energy(a[n:])
    mono={name:energy(read(HERE/'results/clusters'/(name+'.refined.xyz'))) for name in [record['donor'],record['acceptor']]}
    result=dict(pair=args.pair,model_sha256=SHA,geometry=record['best_geometry'],
                frozen_interaction_kJ_mol=(ab-fa-fb)*96.4853321233,
                binding_relative_refined_monomers_kJ_mol=(ab-mono[record['donor']]-mono[record['acceptor']])*96.4853321233,
                scope='Electronic gas-cluster energies, not liquid chemical potentials or cleaning selectivity')
    target.write_text(json.dumps(result,indent=2)+'\n');print(result)


if __name__=='__main__':main()
