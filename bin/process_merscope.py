#!/usr/bin/env python3
"""
process_merscope.py
====================================
Processes MERSCOPE output
"""

from __future__ import annotations
import os
import fire
import logging
import h5py
import scanpy as sc
import numpy as np
import pandas as pd
import tifffile as tf
from skimage.draw import polygon
from process_h5ad import h5ad_to_zarr


def merscope_to_anndata(path: str, filter_prefix: str = "Blank-") -> sc.AnnData:
    """Function to create an AnnData object from a MERSCOPE output directory.

    Args:
        path (str): Path to a MERSCOPE output directory
        filter_prefix (str, optional): Prefix to filter out of the data. Defaults to `Blank-`

    Returns:
        AnnData: AnnData object created from the MERSCOPE output data
    """

    m = pd.read_csv(os.path.join(path, "cell_by_gene.csv"), index_col=0)
    if filter_prefix and len(filter_prefix):
        m = m.loc[:, ~m.columns.str.startswith(filter_prefix)]
    m.index = m.index.astype(str)

    cells = pd.read_csv(os.path.join(path, "cell_metadata.csv"), index_col=0)
    cells.index = cells.index.astype(str)

    adata = sc.AnnData(
        m.to_numpy(),
        dtype=np.float32,
        obs=cells.loc[m.index],
        var=m.columns.to_frame(name="gene"),
    )

    adata.obs["total_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes_by_counts"] = (adata.X > 0).sum(axis=1)

    adata.var["total_counts"] = adata.X.sum(axis=0)
    adata.var["n_cells_by_counts"] = (adata.X > 0).sum(axis=0)

    tm = pd.read_csv(
        os.path.join(path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
        dtype=float,
    ).values

    sp_coords = adata.obs[["center_x", "center_y"]].values
    sp_coords[:, 0] = sp_coords[:, 0] * tm[0, 0] + tm[0, 2]
    sp_coords[:, 1] = sp_coords[:, 1] * tm[1, 1] + tm[1, 2]

    adata.obsm["spatial"] = sp_coords

    return adata


def merscope_to_zarr(
    path: str,
    stem: str,
    filter_prefix: str = "Blank-",
    save_h5ad: bool = False,
    **kwargs,
) -> str:
    """Function to write to Zarr an AnnData object created from MERSCOPE output data

    Args:
        path (str): Path to a MERSCOPE output directory
        stem (str): Prefix for the output Zarr filename
        filter_prefix (str, optional): Prefix to filter out of the data. Defaults to `Blank-`
        save_h5ad (bool, optional): If the AnnData object should also be written to an h5ad file. Defaults to False.

    Returns:
        str: Output Zarr filename
    """

    adata = merscope_to_anndata(path, filter_prefix)
    if save_h5ad:
        adata.write_h5ad(f"{stem}.h5ad")
    zarr_file = h5ad_to_zarr(adata=adata, stem=stem, **kwargs)

    return zarr_file


def merscope_label(
    stem: str, path: str, shape: tuple[int, int], z_index: list[int] = [0]
):
    """This function writes a label image tif file with drawn labels according
    to `cell_boundaries` data stored in MERSCOPE output directory

    Args:
        stem (str): Prefix for the output image filename.
        path (str): Path to the MERSCOPE output directory
        shape (tuple[int, int]): Output image shape.
        z_index (list[int], optional): Z indices to process. Defaults to [0].
    """

    z_index = [z_index] if not isinstance(z_index, (list, tuple)) else z_index

    tm = pd.read_csv(
        os.path.join(path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
        dtype=float,
    ).values

    fovs = [
        x
        for x in os.listdir(os.path.join(path, "cell_boundaries"))
        if x.endswith(".hdf5")
    ]

    for i, z in [(x, "zIndex_{}".format(x)) for x in z_index]:

        label_img = np.zeros((shape[0], shape[1]), dtype=np.uint32)

        for fov in fovs:
            with h5py.File(os.path.join(path, "cell_boundaries", fov)) as f:
                for cell_id in f["featuredata"].keys():
                    pol = f["featuredata"][cell_id][z]["p_0"]["coordinates"][0]

                    pol[:, 0] = pol[:, 0] * tm[0, 0] + tm[0, 2]
                    pol[:, 1] = pol[:, 1] * tm[1, 1] + tm[1, 2]

                    rr, cc = polygon(pol[:, 1], pol[:, 0])
                    label_img[rr - 1, cc - 1] = int(cell_id)

        logging.info(f"Writing label tif image {stem}-label_z{i} ...")
        tf.imwrite(
            f"{stem}-label_z{i}.tif" if len(z_index) > 1 else f"{stem}-label.tif",
            label_img,
        )

    return


if __name__ == "__main__":
    fire.Fire(merscope_to_zarr)
