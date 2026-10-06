"""Exact bounded radial network + tensor product + neighbor aggregation.

Do not create all edge weights or their adjoints at once. Each edge uses
the original radial networks, tensor product, cutoff and density model.
First-derivative inference only, frozen parameters; no LAMMPS ghosts.
"""
import torch


def _edge_values(x, y, feats, cutoff, operation):
    weights=operation.conv_tp_weights(feats)*cutoff
    density=torch.tanh(operation.density_fn(feats)**2)*cutoff
    return operation.conv_tp(x,y,weights),density


class _RadialMessage(torch.autograd.Function):
    @staticmethod
    def forward(ctx,x,y,feats,edges,cutoff,operation,chunk):
        ctx.operation,ctx.chunk=operation,chunk
        ctx.save_for_backward(x,y,feats,edges,cutoff)
        assert all(not p.requires_grad for p in operation.parameters())
        message=density=None
        for start in range(0,len(y),chunk):
            stop=min(start+chunk,len(y));source,target=edges[:,start:stop]
            part,rho=_edge_values(x[source],y[start:stop],feats[start:stop],cutoff[start:stop],operation)
            if message is None:
                message=torch.zeros((len(x),*part.shape[1:]),dtype=part.dtype,device=part.device)
                density=torch.zeros((len(x),*rho.shape[1:]),dtype=rho.dtype,device=rho.device)
            message.index_add_(0,target,part);density.index_add_(0,target,rho)
        return message,density

    @staticmethod
    def backward(ctx,grad_message,grad_density):
        if torch.is_grad_enabled():raise RuntimeError('Radial fusion supports first derivatives only')
        x,y,feats,edges,cutoff=ctx.saved_tensors
        slots=(0,1,2,4);needed=[ctx.needs_input_grad[i] for i in slots]
        saved=(x,y,feats,cutoff)
        out=[torch.zeros_like(v) if need else None for v,need in zip(saved,needed)]
        for start in range(0,len(y),ctx.chunk):
            stop=min(start+ctx.chunk,len(y));source,target=edges[:,start:stop]
            with torch.enable_grad():
                values=[v.detach().requires_grad_(need) for v,need in zip(
                    (x[source],y[start:stop],feats[start:stop],cutoff[start:stop]),needed)]
                outputs=_edge_values(*values,ctx.operation)
                active=[v for v,need in zip(values,needed) if need]
                pairs=[(v,g) for v,g in zip(outputs,(grad_message[target],grad_density[target])) if v.requires_grad]
                grads=torch.autograd.grad([v for v,g in pairs],active,grad_outputs=[g for v,g in pairs],
                                          allow_unused=True,create_graph=False,retain_graph=False)
            cursor=0
            for i,need in enumerate(needed):
                if need:
                    g=grads[cursor];cursor+=1
                    if g is None:continue
                    if i==0:out[i].index_add_(0,source,g)
                    else:out[i][start:stop]=g
        return out[0],out[1],out[2],None,out[3],None,None


class FusedRadialInteraction(torch.nn.Module):
    def __init__(self,inner,chunk):
        super().__init__();self.inner,self.chunk=inner,chunk

    def forward(self,node_attrs,node_feats,edge_attrs,edge_feats,edge_index,cutoff=None,
                lammps_class=None,lammps_natoms=(0,0),first_layer=False):
        assert lammps_class is None and lammps_natoms==(0,0) and cutoff is not None
        op=self.inner
        sc=op.skip_tp(node_feats)
        features=op.linear_up(node_feats);residual=op.linear_res(features)
        # Stock handle_lammps is identity with lammps_class=None.
        source_embedding=op.source_embedding(node_attrs)
        target_embedding=op.target_embedding(node_attrs)
        encoded=torch.cat([edge_feats,source_embedding[edge_index[0]],target_embedding[edge_index[1]]],dim=-1)
        message,density=_RadialMessage.apply(features,edge_attrs,encoded,edge_index,cutoff,op,self.chunk)
        message=op.linear_1(message)/(density*op.beta+op.alpha)
        message=op.linear_2(op.equivariant_nonlin(message+residual))
        return op.reshape(message),sc


def wrap_radial_interactions(model,chunk=256):
    assert chunk>0
    changed=[]
    for i,inner in enumerate(model.interactions):
        assert type(inner).__name__=='RealAgnosticResidualNonLinearInteractionBlock'
        assert not hasattr(inner,'conv_fusion')
        model.interactions[i]=FusedRadialInteraction(inner,chunk)
        changed.append(f'interactions.{i}.radial+product+sum_exact_chunks')
    return changed
