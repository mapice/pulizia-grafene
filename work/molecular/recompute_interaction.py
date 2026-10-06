"""First-derivative inference by exact full interaction recomputation.

Store the interaction inputs, not radial-network and edge-weight outputs
of every layer. This avoids checkpoint's partial TorchScript replay.
Frozen weights, no LAMMPS ghosts, no second derivatives/training supported.
"""
import torch


class _Interaction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, attrs, feats, edge_attrs, edge_feats, edges, cutoff, operation, first):
        ctx.operation, ctx.first = operation, first
        ctx.save_for_backward(attrs, feats, edge_attrs, edge_feats, edges, cutoff)
        assert all(not p.requires_grad for p in operation.parameters())
        return operation(node_attrs=attrs, node_feats=feats, edge_attrs=edge_attrs,
                         edge_feats=edge_feats, edge_index=edges, cutoff=cutoff, first_layer=first)

    @staticmethod
    def backward(ctx, grad_feats, grad_skip):
        if torch.is_grad_enabled():
            raise RuntimeError('Interaction recomputation supports first derivatives only')
        needed = ctx.needs_input_grad[:6]
        with torch.enable_grad():
            values = [v.detach().requires_grad_(need) for v, need in zip(ctx.saved_tensors, needed)]
            outputs = ctx.operation(node_attrs=values[0], node_feats=values[1], edge_attrs=values[2],
                                    edge_feats=values[3], edge_index=values[4], cutoff=values[5], first_layer=ctx.first)
            active = [v for v, need in zip(values, needed) if need]
            pairs = [(v, g) for v, g in zip(outputs, (grad_feats, grad_skip)) if v.requires_grad and g is not None]
            if active and pairs:
                gradients = torch.autograd.grad([v for v, _ in pairs], active,
                                                grad_outputs=[g for _, g in pairs], allow_unused=True,
                                                create_graph=False, retain_graph=False)
            else:
                gradients = [None] * len(active)
        result, cursor = [], 0
        for need in needed:
            result.append(gradients[cursor] if need else None)
            cursor += int(need)
        return *result, None, None


class RecomputeInteraction(torch.nn.Module):
    def __init__(self, inner):
        super().__init__()
        self.inner = inner

    def forward(self, node_attrs, node_feats, edge_attrs, edge_feats, edge_index,
                cutoff=None, lammps_class=None, lammps_natoms=(0, 0), first_layer=False):
        assert lammps_class is None and lammps_natoms == (0, 0) and cutoff is not None
        return _Interaction.apply(node_attrs, node_feats, edge_attrs, edge_feats, edge_index,
                                  cutoff, self.inner, first_layer)


def wrap_interactions(model):
    changed = []
    for i, inner in enumerate(model.interactions):
        assert type(inner).__name__ == 'RealAgnosticResidualNonLinearInteractionBlock'
        model.interactions[i] = RecomputeInteraction(inner)
        changed.append(f'interactions.{i}.exact_recompute_first_derivative')
    return changed
