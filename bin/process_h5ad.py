#!/usr/bin/env python3

import fire
import scanpy as sc
from scipy.sparse import spmatrix
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

SUFFIX = "anndata.zarr"

def h5ad_to_zarr(
    adata,
    stem,
    compute_embeddings=False,
    chunk_size=10,
    ):
    #check if index is integer, if not reindex

    if not adata.obs.index.is_integer(): 
        adata.obs['label_id'] = adata.obs.index
        adata.obs.index = pd.Categorical(adata.obs.index)
        adata.obs.index = adata.obs.index.codes
        adata.obs.index = adata.obs.index.astype(str)
    #turn obsm into a numpy array
    for k in adata.obsm_keys():
        adata.obsm[k] = np.array(adata.obsm[k])
    # compute embeddings if not already stored in object
    if compute_embeddings:
        if not "X_pca" in adata.obsm:
            sc.tl.pca(adata)
        if not "X_umap" in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

    for col in adata.obs:
        # anndata >= 0.8.0
        # if data type is categorical vitessce will throw "path obs/X contains a group" and won"t find .zarray
        # if adata.obs[col].dtype == "category":
        #     adata.obs[col] = adata.obs[col].cat.codes
        if adata.obs[col].dtype == "int8":
            adata.obs[col] = adata.obs[col].astype("int32")
        if adata.obs[col].dtype == "bool":
            adata.obs[col] = adata.obs[col].astype(str).astype("category")

    for col in adata.obsm:
        if type(adata.obsm[col]).__name__ in ["DataFrame", "Series"]:
            adata.obsm[col] = adata.obsm[col].to_numpy()
        if adata.obsm[col].dtype == "int8" or col == "spatial":
            adata.obsm[col] = adata.obsm[col].astype("int32")

    # matrix sparse to dense
    if isinstance(adata.X, spmatrix):
        # use toarray() as it generates a ndarray, instead of todense() which generates a matrix
        adata.X = adata.X.toarray()

    # remove unnecessary data
    del adata.raw

    zarr_file = f"{stem}_{SUFFIX}"
    adata.write_zarr(zarr_file, [adata.shape[0], chunk_size])

    return zarr_file

if __name__ == "__main__":
   fire.Fire(h5ad_to_zarr)
