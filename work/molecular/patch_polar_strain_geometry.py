"""Complete strain derivatives of the unchanged polar-model energy.

The stock graph context strains short-range vectors, while returning the
original positions/cell. Electrostatic positions, reciprocal cell and
volume therefore need their own differentiable strained geometry.
At zero strain energies and position forces must stay invariant; stress
must agree with independent finite cell derivatives in CPU64.
"""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
target=HERE/'.venv/lib/python3.12/site-packages/mace/modules/extensions.py'
root=HERE/'patches';root.mkdir(exist_ok=True)
backup=root/'polar-extensions-pre-strain.py'
text=target.read_text();prefix,body=text.split('class PolarMACE(',1)
marker='''        # Build k-grid
'''
insert='''        # Complete physical strain for the explicit electrostatic geometry.
        # Keep `positions` as the original force differentiation variable.
        field_positions = data["positions"]
        field_cell = cell.view(-1, 3, 3)
        if compute_virials or compute_stress or compute_displacement:
            symmetric_strain = 0.5 * (displacement + displacement.transpose(-1, -2))
            field_cell = field_cell + torch.matmul(field_cell, symmetric_strain)
        # 3x3 linear algebra on CPU for Metal; differentiable device copies,
        # without global fallback settings or double tensors on the GPU.
        geometry_cell = field_cell.to("cpu") if field_cell.device.type == "mps" else field_cell
        geometry_volume = torch.linalg.det(geometry_cell)
        has_cell = geometry_volume.detach().abs() > 1e-12
        identity = torch.eye(3, dtype=geometry_cell.dtype, device=geometry_cell.device)
        safe_cell = torch.where(has_cell[:, None, None], geometry_cell, identity)
        reciprocal = 2 * torch.pi * torch.linalg.inv(safe_cell.transpose(-1, -2))
        reciprocal = torch.where(has_cell[:, None, None], reciprocal, torch.zeros_like(reciprocal))
        field_reciprocal = reciprocal.to(field_cell.device)
        field_volume = geometry_volume.to(field_cell.device)

        # Build k-grid
'''
if 'Complete physical strain for the explicit electrostatic geometry.' in body:
    assert backup.exists()
else:
    assert body.count(marker)==1 and not backup.exists()
    backup.write_text(text)
    body=body.replace(marker,insert)
    old='self.kspace_cutoff, cell.view(-1, 3, 3), data["rcell"].view(-1, 3, 3)'
    assert body.count(old)==1
    body=body.replace(old,'self.kspace_cutoff, field_cell, field_reciprocal')
    assert body.count('node_positions=positions,')==2
    body=body.replace('node_positions=positions,','node_positions=field_positions,')
    assert body.count('volume=data["volume"],')==2
    body=body.replace('volume=data["volume"],','volume=field_volume,')
    body=body.replace('src=_polar_accumulation_dtype(positions),','src=_polar_accumulation_dtype(field_positions),')
    body=body.replace('positions - barycenter[data["batch"], :],','field_positions - barycenter[data["batch"], :],')
    body=body.replace('charge_density_mul_ir, positions, data["batch"], num_graphs',
                      'charge_density_mul_ir, field_positions, data["batch"], num_graphs')
    target.write_text(prefix+'class PolarMACE('+body)
record=dict(before_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),
            after_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
            scope='Derivative consistency correction in isolated environment; neural weights and zero-strain energy preserved',
            independent_verification_required=True)
(root/'polar-strain-geometry.json').write_text(json.dumps(record,indent=2)+'\n')
print(record)
