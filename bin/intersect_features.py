from __future__ import annotations
import typing as T
import os
import h5py
import zarr
import anndata as ad
import pandas as pd
from scipy.sparse import spmatrix


def intersect_features(
    anndata_files: list[str],
    chunk_size: int = 10,
    override_chunk_size: bool = False,
):

    var_indices = []
    for file in anndata_files:
        is_zarr = os.path.splitext(file)[-1] == ".zarr"
        if is_zarr:
            z = zarr.open(file)
            var_indices.append(pd.Index(z.var._index[:]).to_series())
        else:
            with h5py.File(file, "r") as f:
                var_indices.append(ad._io.h5ad.read_elem(f["var"]).index.to_series())

    var_intersect = pd.concat(var_indices, axis=1, join="inner").index

    for file in anndata_files:
        basename, ext = os.path.splitext(file)
        is_zarr = ext == ".zarr"
        if is_zarr:
            z = zarr.open(file)
            adata = ad.read_zarr(z.store)
        else:
            adata = ad.read(file)

        adata = adata[:, var_intersect]

        if is_zarr and not override_chunk_size:
            chunk_shape = z.X.chunks
        else:
            chunk_shape = (adata.shape[0], chunk_size)

        if isinstance(adata.X, spmatrix):
            adata.X = adata.X.toarray()

        adata.write_zarr(f"reindexed-{basename}.zarr", chunk_shape)

    return
