"""Exact first-derivative message passing without storing every edge message.

All original edges and model weights are used. The tensor product and the
sum over target atoms are evaluated together, in bounded chunks. Frozen
model inference only; second derivatives/training are explicitly rejected.
"""
import torch


class _MessageFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, y, w, edges, operation, chunk):
        assert y.shape[0] == w.shape[0] == edges.shape[1]
        ctx.operation, ctx.chunk = operation, chunk
        ctx.save_for_backward(x, y, w, edges)
        out = None
        for start in range(0, len(y), chunk):
            stop = min(start + chunk, len(y))
            part = operation(x[edges[0, start:stop]], y[start:stop], w[start:stop])
            if out is None:
                out = torch.zeros((len(x), *part.shape[1:]), dtype=part.dtype, device=part.device)
            out.index_add_(0, edges[1, start:stop], part)
        assert out is not None
        return out

    @staticmethod
    def backward(ctx, grad_output):
        if torch.is_grad_enabled():
            raise RuntimeError('Fused inference supports only first derivatives, not Hessians or training')
        x, y, w, edges = ctx.saved_tensors
        needed = ctx.needs_input_grad[:3]
        out = [torch.zeros_like(v) if need else None for v, need in zip((x, y, w), needed)]
        for start in range(0, len(y), ctx.chunk):
            stop = min(start + ctx.chunk, len(y))
            source, target = edges[:, start:stop]
            with torch.enable_grad():
                args = [v.detach().requires_grad_(need) for v, need in zip(
                    (x[source], y[start:stop], w[start:stop]), needed)]
                value = ctx.operation(*args)
                active = [v for v, need in zip(args, needed) if need]
                if not active:
                    continue
                gradients = torch.autograd.grad(value, active, grad_outputs=grad_output[target],
                                                allow_unused=True, create_graph=False, retain_graph=False)
            cursor = 0
            for i, need in enumerate(needed):
                if need:
                    gradient = gradients[cursor]
                    cursor += 1
                    if gradient is None:
                        continue
                    if i == 0:
                        out[i].index_add_(0, source, gradient)
                    else:
                        out[i][start:stop] = gradient
        return *out, None, None, None


class FusedChunkMessage(torch.nn.Module):
    def __init__(self, inner, chunk):
        super().__init__()
        self.inner, self.chunk = inner, chunk

    def forward(self, x, y, w, edge_index):
        return _MessageFunction.apply(x, y, w, edge_index, self.inner, self.chunk)


def wrap_messages(model, chunk=256):
    assert chunk > 0
    changed = []
    for i, interaction in enumerate(model.interactions):
        assert not hasattr(interaction, 'conv_fusion')
        assert type(interaction).__name__ == 'RealAgnosticResidualNonLinearInteractionBlock'
        interaction.conv_tp = FusedChunkMessage(interaction.conv_tp, chunk)
        # The stock forward method already contains this exact fused interface.
        interaction.conv_fusion = True
        changed.append(f'interactions.{i}.conv_tp+target_sum')
    return changed
