"""Evaluate official ByteFF-Pol on the same gas proxies, safely on CPU.

Weights-only checkpoint loading; immutable source/model revisions. No vendor
installation shell script, global OpenMM patch or system configuration change.
Passing molecular tests does not qualify graphene, silica or liquid cleaning.
"""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[v] = '1'
os.nice(8)
from pathlib import Path
import argparse, hashlib, json, sys, time
import numpy as np
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float64)
from rdkit import Chem
from ase.io import read

HERE = Path(__file__).resolve().parent
VENDOR = HERE / 'vendor/byteff-pol'
sys.path[:0] = [str(VENDOR), str(VENDOR / 'submodules/bytemol')]
from byteff2.train.utils import load_model, get_nb_params
from byteff2.data import ClusterData, collate_data
from bytemol.core import Molecule
MOLECULES = {name: value['smiles'] for name, value in
             json.loads((HERE / 'results/clusters/monomers.json').read_text()).items()}
MOLECULES['Anisole'] = 'COc1ccccc1'

HASH = 'ae47a6e6860b563908a2e0a83d4a3f6adc1c36b48f544e2241d24066d43d539c'


def mapped_molecule(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    for i, a in enumerate(mol.GetAtoms(), 1):
        a.SetAtomMapNum(i)
    text = Chem.MolToSmiles(mol)
    obj = Molecule.from_mapped_smiles(text, nconfs=0)
    assert np.array_equal(obj.atomic_numbers, [a.GetAtomicNum() for a in mol.GetAtoms()])
    # Map identifiers must preserve every bond, beyond element-order checks.
    remapped = Chem.MolFromSmiles(obj.get_mapped_smiles(), sanitize=False)
    actual = {tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())))
              for b in remapped.GetBonds()}
    expected = {tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum())))
                for b in mol.GetBonds()}
    assert actual == expected
    return obj


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--parameters-only', action='store_true')
    args = p.parse_args()
    source = HERE / 'sources/byteff-pol'
    assert hashlib.sha256((source / 'optimal.pt').read_bytes()).hexdigest() == HASH
    model = load_model(str(source))
    out = HERE / 'results/byteff-pol'
    out.mkdir(exist_ok=True)
    monomers = {}
    start = time.monotonic()
    for name in ['HFIP', 'DMF', 'THF', 'Anisole', 'AceticAcid', 'MethylPivalate', 'DimethylCarbonate']:
        obj = mapped_molecule(MOLECULES[name])
        meta, params, _, _ = get_nb_params(model, obj, relax=False)
        charge = float(sum(params['charge']))
        assert abs(charge) < 1e-7 and np.all(np.asarray(params['alpha']) > 0)
        assert len(params['charge']) == obj.natoms
        record = dict(name=name, mapped_smiles=obj.get_mapped_smiles(),
                      atom_numbers=list(obj.atomic_numbers), atoms=obj.natoms,
                      formal_charge_e=charge, metadata=meta, parameters=params,
                      gas_or_liquid_qualification_not_implied=True,
                      source_commit='8f2813407ba5fbecfb5ec5c69e10b124c5b5bdc2', model_sha256=HASH)
        (out / (name + '-parameters.json')).write_text(json.dumps(record, indent=2) + '\n')
        monomers[name] = obj
        print(name, obj.natoms, 'q', charge, 'sum atomic alpha nm3', sum(params['alpha']), flush=True)
    if args.parameters_only:
        return

    def pair_evaluation(names, coords):
        objects = [monomers.get(n) or mapped_molecule(MOLECULES[n]) for n in names]
        numbers = np.concatenate([np.asarray(m.atomic_numbers) for m in objects])
        data = ClusterData('test', [m.get_mapped_smiles() for m in objects],
                           confdata={'coords': np.asarray(coords)[None, :, :]}, max_n_confs=1,
                           float_dtype=torch.float64)
        data = collate_data([data])
        # Vendor ClusterData constructs coordinates as float32 explicitly.
        # Match all real-valued data to the CPU64 model before evaluation;
        # integer topology/index tensors retain their exact dtype.
        for field, value in list(data.items()):
            if torch.is_tensor(value) and value.is_floating_point():
                data[field] = value.double()
        result = model(data, cluster=True)
        # Same subtraction used by the authors' compare_dimer.py.
        e = float((result['energy_cluster'].sum() - result['energy'].sum()).detach()) * 4.184
        fields = {k: float(result['ff_parameters'][k].sum().detach()) * 4.184
                  for k in ['DISP', 'PAULI', 'ELEC', 'POLARIZATION', 'CHARGE_TRANSFER']}
        assert abs(sum(fields.values()) - e) < 1e-6, (fields, e)
        forces = result['forces_cluster'].detach().numpy().reshape(-1, 3) * 4.184
        assert np.isfinite(e) and np.isfinite(forces).all()
        return e, fields, forces, numbers

    pairs = json.loads((HERE / 'results/clusters/binding-refined.json').read_text())['pairs']
    rows = []
    for key, record in pairs.items():
        names = [record['acceptor'], record['donor']]
        a = read(HERE / 'results/clusters' / record['best_geometry'])
        n = next(t['n_acceptor'] for t in record['trials'] if t['name'] + '.refined.xyz' == record['best_geometry'])
        e, fields, forces, numbers = pair_evaluation(names, a.positions)
        assert np.array_equal(numbers, a.numbers)
        row = dict(pair=key, geometry=record['best_geometry'],
                   interaction_kJ_mol=e, components_kJ_mol=fields,
                   component_sum_independent_error_kJ_mol=sum(fields.values())-e)
        if key == 'HFIP__MethylAcetate':
            h = 1e-4
            plus, minus = a.positions.copy(), a.positions.copy()
            plus[n:, 0] += h
            minus[n:, 0] -= h
            ep, *_ = pair_evaluation(names, plus)
            em, *_ = pair_evaluation(names, minus)
            f_numeric = -(ep-em)/(2*h)
            f_model = float(forces[n:, 0].sum())
            assert abs(f_model-f_numeric) < 2e-4, (f_model, f_numeric)
            far = a.positions.copy()
            far[n:] += [40., 0., 0.]
            efar, *_ = pair_evaluation(names, far)
            assert abs(efar) < .05, efar
            row.update(rigid_fragment_force_finite_difference_kJ_mol_A=f_numeric,
                       rigid_fragment_force_reported_kJ_mol_A=f_model,
                       far_interaction_kJ_mol=efar)
        rows.append(row)
        print(key, e, fields, flush=True)
        (out / 'pairs-partial.json').write_text(json.dumps(rows, indent=2) + '\n')
    result = dict(model='ByteFF-Pol officialv1.0.0', model_sha256=HASH,
                  torch_version=torch.__version__, precision='CPU float64; state_dict loaded with weights_only=True',
                  coordinate_construction_dtype='float64 from numpy input, not cast after float32 rounding',
                  scope='Frozen gas proxies. No extended graphene/silica qualification, no liquid free energy or residue cleaning.',
                  experimental_efficacy_not_established=True, pairs=rows,
                  elapsed_seconds=time.monotonic()-start,
                  source_openmm_patches_not_applied_to_existing_installation=True)
    (out / 'checks.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Finished', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
