#!/usr/bin/env python3
"""
generate_label.py
====================================
Generates the label image from AnnData spatial data
"""

from __future__ import annotations
import os
import fire
import typing as T
import scanpy as sc
import numpy as np
from skimage.draw import disk
import tifffile as tf
import pandas as pd
from pathlib import Path
from process_spaceranger import spaceranger_to_anndata


def from_h5ad(stem: str, h5ad: str, shape: tuple[int, int] = None) -> None:
    sample_id = Path(h5ad).stem

    if os.path.isdir(h5ad):
        adata = spaceranger_to_anndata(h5ad)
        hires_shape = adata.uns["spatial"][sample_id]["images"]["hires"].shape
        scalef = adata.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
        shape = [ hires_shape[0] / scalef, hires_shape[1] / scalef ]
    else:
        adata = sc.read(h5ad)
    
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

    labelImg = np.zeros((shape[0], shape[1]), dtype=np.uint16)

    for spId, (y, x) in zip(adata.obs.index, spot_coords):
        labelImg[disk((x, y), spot_diameter_fullres / 2)] = int(spId)

    tf.imwrite(f"{stem}.tif", labelImg)

    return


def create_img(stem: str, file_type: str, file_path: str, ref_img: str, args: dict[str, T.Any]) -> None:
    """This function writes a label image of spots
    from the metadata contained in an h5ad file.

    Args:
        stem (str): Prefix for the output image filename
        ome_md (dict of str: str): Dictionary containing image size in keys `X` and `Y`
        h5ad (str): Path to h5ad file
    """

    if file_type == "h5ad":
        from_h5ad(stem, file_path, **args)
    

    


if __name__ == "__main__":
    fire.Fire(create_img)
