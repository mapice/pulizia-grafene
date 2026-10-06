"""Compare internal/intermolecular energies after both longer-HMC pilots.

Frozen endpoints, not a stationary ensemble, equation of state or cleaning
prediction. Import the earlier audit to keep its numerical setup identical.
"""
import audit_mc_internal_energy as setup
from pathlib import Path
import hashlib,json,time
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=setup.HERE;ROOT=HERE/'results/hmc-endpoint-energy-audit';KJMOL=setup.KJMOL


def main():
    assert not ROOT.exists();ROOT.mkdir(parents=True);setup.ROOT=ROOT
    base=HERE/'results/neural-molecular-mc'
    names=[('dense','HFIP_N16_seed20261026_1000moves_continue2000_continue200_hmc12-32_p0.050'),
           ('dilute','HFIP_N16_seed20261027_1000moves_continue200_hmc12-32_p0.050')]
    frames=[]
    for label,name in names:
        parent=base/name;assert (parent/'complete.json').exists()
        result=json.loads((parent/'complete.json').read_text())
        assert not result['liquid_equilibrium_qualified']
        a,provenance=setup.snapshot(parent/'frames.xyz',-1,label)
        frames.append((label,a,provenance))
    source=HERE/'sources/MACE-POLAR-1-M.model'
    assert hashlib.sha256(source.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(source,map_location='cpu',weights_only=False).double()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='cpu',default_dtype='float64')

    def energy(a):
        b=a.copy();b.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
        batch=calc._atoms_to_batch(b).to_dict()
        with torch.no_grad():result=model(batch,compute_force=False,compute_stress=False,compute_virials=False,training=False)
        value=float(result['energy'].item())*calc.energy_units_to_eV
        assert np.isfinite(value);return value

    gas=read(HERE/'results/neural-bulk/HFIP_N16_seed20261027_rho1.000/monomer.xyz')
    reference=energy(gas);rows=[];start=time.monotonic()
    for label,a,provenance in frames:
        length=float(a.cell[0,0]);xyz=a.positions.reshape(16,12,3)
        relative=xyz-xyz[:,:1,:];relative-=length*np.rint(relative/length)
        assert np.linalg.norm(relative,axis=2).max()<length/2
        a.positions=(xyz[:,:1,:]+relative).reshape(-1,3)
        total=energy(a);isolated=[]
        for i in range(16):
            b=a[12*i:12*(i+1)];b.pbc=False;b.set_cell(np.zeros((3,3)));b.positions-=b.positions.mean(axis=0)
            isolated.append(energy(b))
        row=dict(branch=label,provenance=provenance,density_g_cm3=float(a.get_masses().sum()/(.602214076*a.get_volume())),
                 total_energy_eV=total,isolated_fragment_energy_eV=isolated,gas_reference_eV=reference,
                 internal_energy_above_gas_kJ_mol_per_molecule=float((np.mean(isolated)-reference)*KJMOL),
                 intermolecular_energy_kJ_mol_per_molecule=float((total-sum(isolated))/16*KJMOL),
                 elapsed_seconds=time.monotonic()-start)
        rows.append(row);(ROOT/'partial.json').write_text(json.dumps(rows,indent=2)+'\n')
        print({k:row[k] for k in ['branch','density_g_cm3','internal_energy_above_gas_kJ_mol_per_molecule','intermolecular_energy_kJ_mol_per_molecule']},flush=True)
    result=dict(model_sha256=SHA,rows=rows,
                internal_energy_branch_difference_kJ_mol_per_molecule=rows[0]['internal_energy_above_gas_kJ_mol_per_molecule']-rows[1]['internal_energy_above_gas_kJ_mol_per_molecule'],
                scope='Two frozen endpoints after short randomized-HMC continuations; no equilibrium or free-energy inference',
                liquid_equilibrium_qualified=False,cleaning_solution_established=False)
    (ROOT/'complete.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
