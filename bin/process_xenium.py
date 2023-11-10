#!/usr/bin/env python3
"""
process_xenium.py
====================================
Processes Xenium output
"""

from __future__ import annotations
import os
import fire
import json
import zarr
import numpy as np
import scanpy as sc
import pandas as pd
import tifffile as tf
from pathlib import Path
from skimage.draw import polygon
from process_h5ad import h5ad_to_zarr


def xenium_to_anndata(
    path: str,
    spatial_as_pixel: bool = True,
    resolution: float = 0.2125,
    load_clusters: bool = True,
    load_embeddings: bool = True,
) -> sc.AnnData:
    """Function to create an AnnData object from Xenium output.

    Args:
        path (str): Path to a xenium output directory
        spatial_as_pixel (bool, optional): Boolean indicating whether spatial coordinates should be
        converted to pixels. Defaults to True.
        resolution (float, optional): Pixel resolution. Defaults to 0.2125.
        load_clusters (bool, optional): If cluster files should be included in the
            AnnData object. Defaults to True.
        load_embeddings (bool, optional): If embedding coordinates files should be included
            in the AnnData object. Defaults to True.

    Returns:
        AnnData: AnnData object created from the xenium output data
    """

    path = Path(path)

    matrix_file = os.path.join(path, "cell_feature_matrix.h5")
    cells_file = os.path.join(path, "cells.csv.gz")

    adata = sc.read_10x_h5(matrix_file)

    with open(os.path.join(path, "experiment.xenium"), "rt") as f:
        adata.uns["xenium"] = json.load(f)

    adata.obs = pd.read_csv(cells_file, compression="gzip", index_col="cell_id")

    adata.obsm["X_spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    if spatial_as_pixel:
        adata.obsm["X_spatial"] = adata.obsm["X_spatial"] / resolution

    if load_clusters:
        for cluster in [
            d for d in (path / "analysis" / "clustering").iterdir() if d.is_dir()
        ]:
            cluster_name = cluster.name.replace("gene_expression_", "")
            cluster_df = pd.read_csv(
                cluster / "clusters.csv",
                index_col="Barcode",
            )

            clusters = cluster_df.reindex(adata.obs.index, fill_value="Undefined")
            adata.obs[cluster_name] = pd.Categorical(clusters["Cluster"].astype(str))

    if load_embeddings:
        embeddings = [
            ("umap", "2_components"),
            ("pca", "10_components"),
        ]
        for embedding, components in embeddings:
            components_name = (
                components
                if (path / "analysis" / embedding / components).exists()
                else f"gene_expression_{components}"
            )
            embedding_df = pd.read_csv(
                os.path.join(
                    path / "analysis" / embedding / components_name / "projection.csv"
                ),
                index_col="Barcode",
            )

            emb = embedding_df.reindex(adata.obs.index, fill_value=0)
            adata.obsm[f"X_{embedding}"] = emb.values

    # starting on v1.3 cell_id looks like "aaabinlp-1"
    # pd.Categorical.codes converts them to int this is done manually at this step
    # instead of reindex_anndata so we control what matches the label image
    adata.obs = adata.obs.reset_index()
    adata.obs.index = (pd.Categorical(adata.obs["cell_id"]).codes + 1).astype(str)

    return adata


def xenium_to_zarr(
    path: str,
    stem: str,
    spatial_as_pixel: bool = True,
    resolution: float = 0.2125,
    save_h5ad: bool = False,
    **kwargs,
) -> str:
    """Function to write to Zarr an AnnData object created from xenium output data

    Args:
        path (str): Path to a xenium output directory
        stem (str): Prefix for the output Zarr filename
        spatial_as_pixel (bool, optional): Boolean indicating whether spatial coordinates should be
        converted to pixels. Defaults to True.
        resolution (float, optional): Pixel resolution. Defaults to 0.2125.
        save_h5ad (bool, optional): If the AnnData object should also be written to an h5ad file. Defaults to False.

    Returns:
        str: Output Zarr filename
    """

    adata = xenium_to_anndata(path, spatial_as_pixel, resolution)
    if save_h5ad:
        adata.write_h5ad(f"tmp-{stem}.h5ad")
    zarr_file = h5ad_to_zarr(adata=adata, stem=stem, **kwargs)

    return zarr_file


def xenium_label(
    stem: str, path: str, shape: tuple[int, int], resolution: float = 0.2125
) -> None:
    """This function writes a label image tif file with drawn labels according to
    cell segmentation polygons from Xenium output cells.zarr.zip file

    Args:
        stem (str): Prefix for the output image filename.
        path (str): Path to the Xenium output directory or cells.zarr.zip file
        shape (tuple[int, int]): Output image shape. Defaults to None.
        resolution (float, optional): Pixel resolution. Defaults to 0.2125.
    """
    if os.path.isdir(path):
        cells_file = os.path.join(path, "cells.zarr.zip")
    else:
        cells_file = path

    z = zarr.open(cells_file, "r")

    with open(os.path.join(path, "experiment.xenium")) as f:
        experiment = json.load(f)
    sw_version = float(experiment["analysis_sw_version"][7:10])

    if sw_version < 1.3:
        ids = z["cell_id"]
    else:
        ids = z["cell_id"][:, 0]

    # starting on v1.3 cell_id looks like "aaabinlp-1"
    # pd.Categorical.codes converts them to int
    # this is required so the label image matches the h5ad ids
    ids = pd.Categorical(ids).codes + 1

    pols = z["polygon_vertices"][1]

    label_img = np.zeros((shape[0], shape[1]), dtype=np.min_scalar_type(max(ids)))

    for id, pol in zip(ids, pols):
        pol = pol / resolution
        pol = np.array(list(map(list, pol.reshape(pol.shape[0] // 2, 2))))
        rr, cc = polygon(pol[:, 1], pol[:, 0])
        label_img[rr - 1, cc - 1] = int(id)

    tf.imwrite(f"{stem}-label.tif", label_img)

    return


if __name__ == "__main__":
    fire.Fire(xenium_to_zarr)
