"""Enable exact reciprocal-basis derivatives required for cell pressure.

Replace out= multiply buffers by the same expressions stacked, and avoid
0/0 from sqrt at the masked zero reciprocal vector. Scalar values stay
unchanged; all tests still require independent CPU64 cell derivatives.
"""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
target=HERE/'.venv/lib/python3.12/site-packages/graph_longrange/gto_utils.py'
root=HERE/'patches';root.mkdir(exist_ok=True)
backup=root/'gto-utils-pre-cell-derivatives.py'
text=target.read_text()
old='''        out = torch.empty(
            (*k_mods.shape, self.num_sigma, 2),
            dtype=k_mods.dtype,
            device=k_mods.device,
        )
        torch.mul(self.pref0, exp_term, out=out[..., 0])
        torch.mul(self.pref1, k_mods.unsqueeze(-1) * exp_term, out=out[..., 1])
        return out'''
new='''        return torch.stack(
            (self.pref0 * exp_term, self.pref1 * k_mods.unsqueeze(-1) * exp_term),
            dim=-1,
        )'''
old_zero='''        k_moduli = torch.sqrt(torch.clamp_min(k_norm2, 0.0))
        k_moduli.masked_fill_(k0_mask > 0.0, 0.0)
        return k_moduli'''
new_zero='''        zero = k0_mask > 0.0
        safe_norm2 = torch.where(zero, torch.ones_like(k_norm2), torch.clamp_min(k_norm2, 0.0))
        k_moduli = torch.sqrt(safe_norm2)
        return torch.where(zero, torch.zeros_like(k_moduli), k_moduli)'''
assert text.count(old)==1 and text.count(old_zero)==1 and not backup.exists()
backup.write_text(text)
target.write_text(text.replace(old,new).replace(old_zero,new_zero))
record=dict(before_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),
            after_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
            scope='Differentiability of same reciprocal Gaussian basis; model weights unchanged; validation pending')
(root/'gto-cell-derivatives.json').write_text(json.dumps(record,indent=2)+'\n')
print(record)
