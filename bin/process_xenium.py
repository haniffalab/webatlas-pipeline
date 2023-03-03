#!/usr/bin/env python3
"""
process_xenium.py
====================================
Processes Xenium output
"""

from __future__ import annotations
import os
import fire
import glob
import zarr
import numpy as np
import scanpy as sc
import pandas as pd
import tifffile as tf
from skimage.draw import polygon
from process_h5ad import h5ad_to_zarr


def xenium_to_anndata(
    path: str, spatial_as_pixel: bool = True, resolution: float = 0.2125
) -> sc.AnnData:
    """Function to create an AnnData object from Xenium output.

    Args:
        path (str): Path to a xenium output directory
        spatial_as_pixel (bool, optional): Boolean indicating whether spatial coordinates should be
        converted to pixels. Defaults to True.
        resolution (float, optional): Pixel resolution. Defaults to 0.2125.

    Returns:
        AnnData: AnnData object created from the xenium output data
    """

    matrix_file = glob.glob(os.path.join(path, "*cell_feature_matrix.h5"))[0]
    cells_file = glob.glob(os.path.join(path, "*cells.csv.gz"))[0]

    adata = sc.read_10x_h5(matrix_file)

    df = pd.read_csv(cells_file, compression="gzip", index_col="cell_id")
    df.index = df.index.astype(str)

    adata.obs = df

    adata.obsm["X_spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    if spatial_as_pixel:
        adata.obsm["X_spatial"] = adata.obsm["X_spatial"] / resolution

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
    stem: str, path: str, shape: tuple[int, int], resolution: int = 0.2125
):
    """This function writes a label image tif file with drawn labels according to
    cell segmentation polygons from Xenium output cells.zarr.zip file

    Args:
        stem (str): Prefix for the output image filename.
        path (str): Path to the Xenium output directory or cells.zarr.zip file
        shape (tuple[int, int]): Output image shape. Defaults to None.
        resolution (int, optional): Pixel resolution. Defaults to 0.2125.
    """
    if os.path.isdir(path):
        cells_file = glob.glob(os.path.join(path, "*cells.zarr.zip"))[0]
    else:
        cells_file = path

    z = zarr.open(cells_file, "r")
    ids = z["cell_id"]
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
