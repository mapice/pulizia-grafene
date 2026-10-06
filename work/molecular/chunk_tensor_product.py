"""Exact chunked first derivatives of frozen e3nn tensor products.

Save only inputs; recompute product chunks during the first backward pass.
This is NOT a smaller/truncated model: every edge and original weight stays.
Only energies, forces and stress are supported; not training or Hessians.
Validation against the unmodified CPU64 calculator is mandatory.
"""
import torch


class _ChunkFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, y, w, operation, chunk):
        assert x.shape[0]==y.shape[0]==w.shape[0]
        ctx.operation,ctx.chunk=operation,chunk
        ctx.save_for_backward(x,y,w)
        out=None
        for start in range(0,len(x),chunk):
            stop=min(start+chunk,len(x))
            part=operation(x[start:stop],y[start:stop],w[start:stop])
            if out is None:
                out=torch.empty((len(x),*part.shape[1:]),dtype=part.dtype,device=part.device)
            out[start:stop]=part
        assert out is not None
        return out

    @staticmethod
    def backward(ctx, grad_output):
        if torch.is_grad_enabled():
            raise RuntimeError('Chunked tensor-product path only supports first derivatives; no training/Hessian')
        saved=ctx.saved_tensors
        needed=ctx.needs_input_grad[:3]
        output=[torch.zeros_like(x) if need else None for x,need in zip(saved,needed)]
        for start in range(0,len(saved[0]),ctx.chunk):
            stop=min(start+ctx.chunk,len(saved[0]))
            with torch.enable_grad():
                args=[x[start:stop].detach().requires_grad_(need) for x,need in zip(saved,needed)]
                value=ctx.operation(*args)
                active=[x for x,need in zip(args,needed) if need]
                if not active:continue
                gradients=torch.autograd.grad(value,active,grad_outputs=grad_output[start:stop],
                                              create_graph=False,retain_graph=False,allow_unused=True)
            cursor=0
            for i,need in enumerate(needed):
                if need:
                    gradient=gradients[cursor];cursor+=1
                    if gradient is not None:output[i][start:stop]=gradient
        return *output,None,None


class ChunkedTensorProduct(torch.nn.Module):
    def __init__(self, inner, chunk=256):
        super().__init__();self.inner=inner;self.chunk=chunk

    def forward(self,x,y,w):
        return _ChunkFunction.apply(x,y,w,self.inner,self.chunk)


def wrap_products(model,chunk=256):
    changed=[]
    for i,interaction in enumerate(model.interactions):
        assert not hasattr(interaction,'conv_fusion'), 'Fused backend needs a different exact interface'
        assert not isinstance(interaction.conv_tp,ChunkedTensorProduct)
        interaction.conv_tp=ChunkedTensorProduct(interaction.conv_tp,chunk)
        changed.append(f'interactions.{i}.conv_tp')
    return changed
