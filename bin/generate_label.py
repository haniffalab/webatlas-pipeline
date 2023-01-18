#!/usr/bin/env python3
"""
generate_label.py
====================================
Generates the label image from AnnData spatial data
"""

from __future__ import annotations
import os
import fire
import scanpy as sc
import numpy as np
from skimage.draw import disk
import tifffile as tf
import pandas as pd
from pathlib import Path
from process_spaceranger import spaceranger_to_anndata


def main(stem: str, ome_md: dict[str, str], h5ad: str) -> None:
    """This function writes a label image of spots
    from the metadata contained in an h5ad file.

    Args:
        stem (str): Prefix for the output image filename
        ome_md (dict of str: str): Dictionary containing image size in keys `X` and `Y`
        h5ad (str): Path to h5ad file
    """
    sample_id = Path(h5ad).stem

    if os.path.isdir(h5ad):
        adata = spaceranger_to_anndata(h5ad)
    else:
        adata = sc.read(h5ad)

    # if multipel sections per slide -------------

    # sample_ids = np.unique(adata.obs[["sample"]])
    # print(sample_ids)
    # md = {}
    # for i in adata.uns["spatial"]:
    # sample_id = re.search(".*_count_\d+_(.+)_mm.*", i).group(1)
    # md[sample_id] = adata.uns["spatial"][i]["scalefactors"]
    # adata.obs["Y"] = adata.obsm["spatial"][:, 0]
    # adata.obs["X"] = adata.obsm["spatial"][:, 1]
    # ----------------------------------------------

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

    X = ome_md["X"]
    Y = ome_md["Y"]
    labelImg = np.zeros((int(Y), int(X)), dtype=np.uint16)

    for spId, (y, x) in zip(adata.obs.index, spot_coords):
        labelImg[disk((x, y), spot_diameter_fullres / 2)] = int(spId)

    tf.imwrite(f"{stem}.tif", labelImg)


if __name__ == "__main__":
    fire.Fire(main)
