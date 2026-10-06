"""Exact algebraic contraction without the batch x channel x channel array.

sum_vd x1[b,u,d]*x2[b,v,d]*w[u,v] is evaluated by first contracting
x2 with w. No approximation, dropped channels or altered weights.
The isolated installed module alone is changed; original source preserved.
Numerical force/stress parity tests against saved CPU64 results required.
"""
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent
target=HERE/'.venv/lib/python3.12/site-packages/mace/modules/field_blocks.py'
root=HERE/'patches';root.mkdir(exist_ok=True)
backup=root/'field-blocks-original.py'
old='''                pair = torch.einsum("bud,bvd->buv", x1v, x2v) / math.sqrt(float(d1))
                mixed = torch.einsum("buv,uv->bu", pair, w)'''
new='''                # Exact reordered contraction avoids B*mul1*mul2 storage.
                weighted_x2 = torch.einsum("bvd,uv->bud", x2v, w)
                mixed = torch.einsum("bud,bud->bu", x1v, weighted_x2) / math.sqrt(float(d1))'''
text=target.read_text()
if new in text:
    assert backup.exists()
else:
    assert text.count(old)==1
    assert not backup.exists()
    backup.write_text(text)
    target.write_text(text.replace(old,new))
record=dict(original_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),
            patched_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
            scope='Exact contraction-order change in isolated environment; model weights unchanged; parity required')
(root/'sparse-contraction.json').write_text(json.dumps(record,indent=2)+'\n')
print(record)
