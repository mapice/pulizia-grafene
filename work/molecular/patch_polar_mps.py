"""Isolated conditional Metal precision patch; CPU equations unchanged.

Changes explicit float64 accumulation only on MPS, which lacks float64.
No model weights or CPU path changed. Numerical force/stress/charge gates
must pass before any trajectory from the patched path is accepted.
"""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
p=HERE/'.venv/lib/python3.12/site-packages/mace/modules/extensions.py'
original=p.read_text();marker='def _polar_accumulation_dtype('
if marker in original:raise RuntimeError('Already patched; preserve current source')
split=original.index('@compile_mode("script")\nclass PolarMACE(')
head,body=original[:split],original[split:]
mapping={
    'fukui_sources.double()':('_polar_accumulation_dtype(fukui_sources)',1),
    'spin_charge_density[:, :, 0].double()':('_polar_accumulation_dtype(spin_charge_density[:, :, 0])',2),
    'positions.double()':('_polar_accumulation_dtype(positions)',1),
    'current_fukui_sources.double()':('_polar_accumulation_dtype(current_fukui_sources)',1),
    'node_e0.clone().double()':('_polar_accumulation_dtype(node_e0.clone())',1),
    'node_inter_es.clone().double()':('_polar_accumulation_dtype(node_inter_es.clone())',1),
}
counts={}
# Replace longer identifiers first; "fukui_sources" is a suffix of
# "current_fukui_sources" and an unconstrained text substitution is wrong.
for old,(new,count) in sorted(mapping.items(),key=lambda x:-len(x[0])):
    actual=body.count(old)
    assert actual==count,(old,actual,count)
    counts[old]=actual;body=body.replace(old,new)
helper='''def _polar_accumulation_dtype(value: torch.Tensor) -> torch.Tensor:
    # Local compatibility: Metal32; exact original double path on CPU.
    return value if value.device.type == "mps" else value.double()


'''
patched=head+helper+body
assert patched.count(marker)==1
back=HERE/'patches/polar-extensions-original.py'
if back.exists():assert back.read_text()==original
else:back.write_text(original)
p.write_text(patched)
record=dict(path=str(p.relative_to(HERE)),original_sha256=hashlib.sha256(original.encode()).hexdigest(),
            patched_sha256=hashlib.sha256(patched.encode()).hexdigest(),replacements=counts,
            model_weights_unchanged=True,CPU_paths_keep_original_double_accumulation=True,
            scope='MPS32 compatibility only, not physics/model qualification',parity_qualification_pending=True)
(HERE/'patches/polar-metal-compatibility.json').write_text(json.dumps(record,indent=2)+'\n')
print(record)
