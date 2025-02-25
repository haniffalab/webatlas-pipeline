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
from process_h5ad import h5ad_to_zarr, subset_anndata


def add_csv_to_adata_obs(adata_obs, annotations_path, annotations_column_index):
    annot_df = pd.read_csv(annotations_path)
    if annotations_column_index in annot_df.columns:
        annot_df.set_index(annotations_column_index, inplace=True)
    else:
        #using simply first column if none were specified
        annot_df.set_index(annot_df.columns[0], inplace=True)
    #adata_obs.set_index('cell_id', inplace=True)
    adata_obs = pd.merge(adata_obs, annot_df, left_index=True, right_index=True, how='left')
    return adata_obs
    


def xenium_to_anndata(
    path: str,
    spatial_as_pixel: bool = True,
    resolution: float = 0.2125,
    load_clusters: bool = True,
    load_embeddings: bool = True,
    annotations: str = None,
    annotations_column_index: str = None,
    obs_subset: tuple[int, T.Any] = None,
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

    if annotations:
        adata.obs = add_csv_to_adata_obs(adata.obs, annotations, annotations_column_index)
    
    # starting on v1.3 cell_id looks like "aaabinlp-1"
    # pd.Categorical.codes converts them to int this is done manually at this step
    # instead of reindex_anndata so we control what matches the label image
    #just to make sure the column with name as index name column does not exist in adata.obs
    index_name = adata.obs.index.name
    if index_name in adata.obs.columns:
        # Drop the existing column with the same name
        adata.obs = adata.obs.drop(columns=[index_name])
    
    adata.obs = adata.obs.reset_index()
    adata.obs.index = (pd.Categorical(adata.obs.index).codes + 1).astype(str)
        

    
    #important to do subsetting after new reindexing
    if obs_subset:
        adata = subset_anndata(adata = adata, obs_subset = obs_subset)
    
    return adata
    
    

def xenium_to_zarr(
    path: str,
    stem: str,
    spatial_as_pixel: bool = True,
    resolution: float = 0.2125,
    save_h5ad: bool = False,
    annotations: str = None,
    annotations_column_index: str = None,
    obs_subset: tuple[int, T.Any] = None,
    **kwargs,
) -> str:
    """Function to write to Zarr an AnnData object created from xenium output data

    Args:
        path (str): Path to a xenium output directory[
        stem (str): Prefix for the output Zarr filename
        spatial_as_pixel (bool, optional): Boolean indicating whether spatial coordinates should be
        converted to pixels. Defaults to True.
        resolution (float, optional): Pixel resolution. Defaults to 0.2125.
        save_h5ad (bool, optional): If the AnnData object should also be written to an h5ad file. Defaults to False.
        annotations (str): path to csv file with annotations or any additional information for "obs"
        annotations_column_index (str): name of column in the csv with barcodes information of visium spots
        obs_subset (tuple(str, T.Any), optional): Tuple containing an `obs` column name and one or more values
            to use to subset the AnnData object. Defaults to None.
    Returns:
        str: Output Zarr filename
    """
    
    adata = xenium_to_anndata(path, spatial_as_pixel = spatial_as_pixel, resolution = resolution, annotations = annotations, annotations_column_index = annotations_column_index, obs_subset = obs_subset)
    if save_h5ad:
        adata.write_h5ad(f"tmp-{stem}.h5ad")
    
    zarr_file = h5ad_to_zarr(adata=adata, stem=stem, **kwargs)

    return zarr_file


def xenium_label(
    stem: str, 
    path: str, 
    shape: tuple[int, int], 
    resolution: float = 0.2125, 
    obs_subset: tuple[str, T.Any] = None,
    var_subset: tuple[str, T.Any] = None,
    annotations: str = None,
    annotations_column_index: str = None,
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
    
    if obs_subset:
        #firstly prepare dict of cell numbers (new) correponding to (old)
        cells_file = os.path.join(path, "cells.csv.gz")
        adata_obs = pd.read_csv(cells_file, compression="gzip", index_col="cell_id")
        #if annotations != 'None':

        if annotations:
            adata_obs = add_csv_to_adata_obs(adata_obs, annotations, annotations_column_index)
        
        #reindexing (as in original anndata)
        adata_obs = adata_obs.reset_index()
        adata_obs.index = (pd.Categorical(adata_obs["cell_id"]).codes + 1).astype(str)
        adata_obs_sub = adata_obs[adata_obs[obs_subset[0]].isin(obs_subset[1])]
        ids = adata_obs_sub.index.to_numpy()
    else:
        if sw_version < 1.3:
            ids = z["cell_id"]
        else:
            ids = z["cell_id"][:, 0]
        ids = pd.Categorical(ids).codes + 1

    
    # starting on v1.3 cell_id looks like "aaabinlp-1"
    # pd.Categorical.codes converts them to int
    # this is required so the label image matches the h5ad ids
    #ids = pd.Categorical(ids).codes + 1
    #we keep the same ids as in original dataset
    ids = np.array(ids).astype('uint32')
    # starting on v2.0 vertices change location
    if sw_version < 2.0:
        pols = z["polygon_vertices"][1]
    else:
        pols = z["polygon_sets"][1]["vertices"]
    label_img = np.zeros((shape[0], shape[1]), dtype=np.min_scalar_type(max(ids)))
    for idd in ids:
        pol = pols[idd-1]
        pol = pol / resolution
        pol = np.array(list(map(list, pol.reshape(pol.shape[0] // 2, 2))))
        rr, cc = polygon(pol[:, 1], pol[:, 0])
        label_img[rr - 1, cc - 1] = int(idd)
    for id, pol in zip(ids, pols):
        pol = pol / resolution
        pol = np.array(list(map(list, pol.reshape(pol.shape[0] // 2, 2))))
        rr, cc = polygon(pol[:, 1], pol[:, 0])
        label_img[rr - 1, cc - 1] = int(id)
    tf.imwrite(f"{stem}-label.tif", label_img)

    return


if __name__ == "__main__":
    fire.Fire(xenium_to_zarr)
