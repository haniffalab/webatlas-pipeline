#!/usr/bin/env python3

import fire
import os

import anndata as ad
import scanpy as sc
import scipy

import warnings
warnings.filterwarnings("ignore")

SUFFIX = 'anndata.zarr'

def h5ad_to_zarr(
    file,
    stem,
    compute_embeddings=False,
    chunk_size=10,
    ):

    adata = ad.read(file)

    # compute embeddings if not already stored in object
    if compute_embeddings:
        if not 'X_pca' in adata.obsm:
            sc.tl.pca(adata)
        if not 'X_umap' in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

    for col in adata.obs:
        # anndata >= 0.8.0
        # if data type is categorical vitessce will throw "path obs/X contains a group" and won't find .zarray
        # if adata.obs[col].dtype == 'category':
        #     adata.obs[col] = adata.obs[col].cat.codes
        if adata.obs[col].dtype == 'int8':
            adata.obs[col] = adata.obs[col].astype('int32')

    for col in adata.obsm:
        if col == 'spatial' or adata.obsm[col].dtype == 'int8':
            adata.obsm[col] = adata.obsm[col].astype('int32')

    # matrix sparse to dense
    if isinstance(adata.X, scipy.sparse.spmatrix):
        # use toarray() as it generates a ndarray, instead of todense() which generates a matrix
        adata.X = adata.X.toarray()

    # remove unnecessary data
    del adata.raw

    zarr_file = f"{stem}_{SUFFIX}"
    adata.write_zarr(zarr_file, [adata.shape[0], chunk_size])

    return zarr_file

if __name__ == "__main__":
   fire.Fire(h5ad_to_zarr)
