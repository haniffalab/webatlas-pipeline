#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import scanpy as sc
import ome_zarr
import numpy as np
import re
from skimage.draw import disk
import tifffile as tf

from pathlib import Path


def main(stem, ome_md, h5ad):
    sample_id = Path(h5ad).stem

    adata = sc.read(h5ad)
    # print(adata)

    # if multipel sections per slide -------------

    # sample_ids = np.unique(adata.obs[["sample"]])
    # print(sample_ids)
    # md = {}
    # for i in adata.uns["spatial"]:
        # sample_id = re.search(".*_count_\d+_(.+)_mm.*", i).group(1)
        # md[sample_id] = adata.uns["spatial"][i]["scalefactors"]
    # adata.obs["Y"] = adata.obsm["spatial"][:, 0]
    # adata.obs["X"] = adata.obsm["spatial"][:, 1]
    #----------------------------------------------

    # print(adata.uns["spatial"][sample_id]["scalefactors"])

    spot_diameter_fullres = adata.uns["spatial"][sample_id]["scalefactors"]["spot_diameter_fullres"]
    # hires_scalef = adata.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
    spot_coords = adata.obsm["spatial"]
    # print(spot_coords, X, Y, Z, C, T)
    assert adata.obs.shape[0] == spot_coords.shape[0]
    label_ids = np.arange(1, adata.obs.shape[0] + 1)
    label_id_col_name = "label_id"
    if label_id_col_name not in adata.obs.columns:
        adata.obs[label_id_col_name] = label_ids

    X = ome_md["X"]
    Y = ome_md["Y"]
    labelImg = np.zeros((int(Y), int(X)), dtype=np.uint16)
    # print(labelImg.shape)
    # print(np.max(spot_coords))
    for spId, (y, x) in zip(adata.obs[label_id_col_name], spot_coords):
        labelImg[disk((y, x), spot_diameter_fullres/2)] = spId

    tf.imwrite(f"{stem}.tif", labelImg)


if __name__ == "__main__":
    fire.Fire(main)
