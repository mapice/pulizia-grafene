"""Recompute inference activations to reduce autograd memory, same equations."""
import torch
from torch.utils.checkpoint import checkpoint


class CheckpointModule(torch.nn.Module):
    def __init__(self, inner):
        super().__init__();self.inner=inner

    def forward(self, *args, **kwargs):
        return checkpoint(self.inner, *args, use_reentrant=False, **kwargs)


def wrap_model(model):
    changed=[]
    for name in ['interactions','products','field_dependent_charges_maps']:
        modules=getattr(model,name)
        for i, module in enumerate(modules):
            assert not isinstance(module,CheckpointModule)
            modules[i]=CheckpointModule(module);changed.append(f'{name}.{i}')
    model.local_electron_energy=CheckpointModule(model.local_electron_energy)
    changed.append('local_electron_energy')
    return changed
