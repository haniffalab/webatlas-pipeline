#!/usr/bin/env python3

import os
import fire
import scanpy as sc
import pandas as pd
from process_h5ad import h5ad_to_zarr


def spaceranger_to_h5ad(
    path, load_clusters=False, load_embeddings=False, clustering="graphclust"
):

    adata = sc.read_visium(path)

    if load_clusters:
        cluster_df = pd.read_csv(
            os.path.join(path, f"analysis/clustering/{clustering}/clusters.csv"),
            index_col="Barcode",
        )
        adata.obs[f"{clustering}"] = cluster_df.Cluster.astype("category")

    if load_embeddings:
        umap = pd.read_csv(
            os.path.join(path, "analysis/umap/2_components/projection.csv"),
            index_col="Barcode",
        )
        tsne = pd.read_csv(
            os.path.join(path, "analysis/tsne/2_components/projection.csv"),
            index_col="Barcode",
        )
        pca = pd.read_csv(
            os.path.join(path, "analysis/pca/10_components/projection.csv"),
            index_col="Barcode",
        )
        adata.obsm["X_umap"] = umap.values
        adata.obsm["X_tsne"] = tsne.values
        adata.obsm["X_pca"] = pca.values

    return adata


def spaceranger_to_zarr(
    path, stem, load_clusters=False, load_embeddings=False, save_h5ad=False, **kwargs
):

    adata = spaceranger_to_h5ad(path, load_clusters, load_embeddings)
    if save_h5ad:
        adata.write_h5ad(f"{stem}.h5ad")
    zarr_file = h5ad_to_zarr(adata=adata, stem=stem, **kwargs)

    return zarr_file


if __name__ == "__main__":
    fire.Fire(spaceranger_to_zarr)
