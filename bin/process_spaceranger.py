#!/usr/bin/env python3

import os
import fire
import scanpy as sc
import pandas as pd
from process_h5ad import h5ad_to_zarr


def spaceranger_to_anndata(
    path: str,
    load_clusters: bool = False,
    load_embeddings: bool = False,
    clustering: str = "graphclust",
) -> sc.AnnData:
    """Function to create an AnnData object from a SpaceRanger output directory.

    Args:
        path (str): Path to a SpaceRanger output directory
        load_clusters (bool, optional): If cluster files should be included in the
            AnnData object. Defaults to False.
        load_embeddings (bool, optional): If embedding coordinates files should be included
            in the AnnData object. Defaults to False.
        clustering (str, optional): The clustering algorithm to include in the AnnData
            object if `load_clusters` is True. Defaults to "graphclust".

    Returns:
        AnnData: AnnData object created from the SpaceRanger output data
    """

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
    path: str,
    stem: str,
    load_clusters: bool = True,
    load_embeddings: bool = True,
    save_h5ad: bool = False,
    **kwargs,
) -> str:
    """Function to write to Zarr an AnnData object created from SpaceRanger output data

    Args:
        path (str): Path to a SpaceRanger output directory
        stem (str): Prefix for the output Zarr filename
        load_clusters (bool, optional): If cluster files should be included in the
            AnnData object. Defaults to False.
        load_embeddings (bool, optional): If embedding coordinates files should be included
            in the AnnData object. Defaults to False.
        save_h5ad (bool, optional): If the AnnData object should also be written to an h5ad file. Defaults to False.

    Returns:
        str: Output Zarr filename
    """

    adata = spaceranger_to_anndata(path, load_clusters, load_embeddings)
    if save_h5ad:
        adata.write_h5ad(f"{stem}.h5ad")
    zarr_file = h5ad_to_zarr(adata=adata, stem=stem, **kwargs)

    return zarr_file


if __name__ == "__main__":
    fire.Fire(spaceranger_to_zarr)
