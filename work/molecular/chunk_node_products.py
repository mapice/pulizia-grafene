"""Chunk only products that are pointwise in the atom index, same equations."""
import torch
from chunk_tensor_product import _ChunkFunction


class ChunkNodeProduct(torch.nn.Module):
    def __init__(self, inner, chunk=16):
        super().__init__()
        self.inner, self.chunk = inner, chunk

    def forward(self, node_feats, sc, node_attrs):
        assert sc is not None and len(node_feats) == len(sc) == len(node_attrs)
        return _ChunkFunction.apply(node_feats, sc, node_attrs, self.inner, self.chunk)


def wrap_node_products(model, chunk=16):
    assert chunk > 0
    changed = []
    for i, product in enumerate(model.products):
        assert type(product).__name__ == 'EquivariantProductBasisBlock'
        # Source audit: symmetric contractions and linear/skip operations
        # act separately at each atom, with no reduction over atom indices.
        model.products[i] = ChunkNodeProduct(product, chunk)
        changed.append(f'products.{i}.pointwise_atom_chunks')
    return changed
