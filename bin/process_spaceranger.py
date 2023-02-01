#!/usr/bin/env python3
"""
process_spaceranger.py
====================================
Processes SpaceRanger output
"""

from __future__ import annotations
import os
import fire
import typing as T
import numpy as np
import scanpy as sc
import pandas as pd
import tifffile as tf
from skimage.draw import disk
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
        embeddings = [
            ("umap", "2_components"),
            ("tsne", "2_components"),
            ("pca", "10_components"),
        ]
        for embedding, components in embeddings:
            emb = pd.read_csv(
                os.path.join(path, f"analysis/{embedding}/{components}/projection.csv"),
                index_col="Barcode",
            )
            adata.obsm[f"X_{embedding}"] = emb.values

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


def visium_label(
    stem: str,
    file_path: str,
    shape: tuple[int, int] = None,
    obs_subset: tuple[int, T.Any] = None,
    sample_id: str = None,
) -> None:
    """This function writes a label image tif file with drawn labels according to an
    Anndata object with necessary metadata stored within `uns["spatial"]`.

    Args:
        stem (str): Prefix for the output image filename.
        file_path (str): Path to the h5ad file or spaceranger output directory.
        shape (tuple[int, int], optional): Output image shape. Defaults to None.
        obs_subset (tuple(str, T.Any), optional): Tuple containing an `obs` column name and one or more values
            to use to subset the AnnData object. Defaults to None.
        sample_id (str, optional): Sample ID string within the Anndata object. Defaults to None.
    """
    # sample_id = sample_id or Path(file_path).stem

    if os.path.isdir(file_path):
        adata = spaceranger_to_anndata(file_path)
        sample_id = sample_id or list(adata.uns["spatial"].keys())[0]
        if not shape:
            hires_shape = adata.uns["spatial"][sample_id]["images"]["hires"].shape
            scalef = adata.uns["spatial"][sample_id]["scalefactors"][
                "tissue_hires_scalef"
            ]
            shape = [int(hires_shape[0] / scalef), int(hires_shape[1] / scalef)]
    else:
        adata = sc.read(file_path)
        sample_id = sample_id or list(adata.uns["spatial"].keys())[0]

    # Subset adata by obs
    if obs_subset:
        obs_subset[1] = (
            [obs_subset[1]]
            if not isinstance(obs_subset[1], (list, tuple))
            else obs_subset[1]
        )
        adata = adata[adata.obs[obs_subset[0]].isin(obs_subset[1])]

    # check if index is numerical, if not reindex
    if not adata.obs.index.is_integer() and not (
        adata.obs.index.is_object() and all(adata.obs.index.str.isnumeric())
    ):
        adata.obs["label_id"] = adata.obs.index
        adata.obs.index = pd.Categorical(adata.obs.index)
        adata.obs.index = adata.obs.index.codes
        adata.obs.index = adata.obs.index.astype(str)

    # turn obsm into a numpy array
    for k in adata.obsm_keys():
        adata.obsm[k] = np.array(adata.obsm[k])

    spot_diameter_fullres = adata.uns["spatial"][sample_id]["scalefactors"][
        "spot_diameter_fullres"
    ]
    # hires_scalef = adata.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
    spot_coords = adata.obsm["spatial"]
    assert adata.obs.shape[0] == spot_coords.shape[0]

    label_img = np.zeros((shape[0], shape[1]), dtype=np.uint16)

    for spId, (y, x) in zip(adata.obs.index, spot_coords):
        label_img[disk((x, y), spot_diameter_fullres / 2)] = int(spId)

    tf.imwrite(f"{stem}-label.tif", label_img)

    return


if __name__ == "__main__":
    fire.Fire(spaceranger_to_zarr)
