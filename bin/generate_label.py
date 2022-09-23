#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import os
import fire
import scanpy as sc
import numpy as np
from skimage.draw import disk
import tifffile as tf
import pandas as pd
from pathlib import Path
from process_spaceranger import spaceranger_to_h5ad


def main(stem, ome_md, h5ad):
    sample_id = Path(h5ad).stem

    if os.path.isdir(h5ad):
        adata = spaceranger_to_h5ad(h5ad)
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

    # check if index is integer, if not reindex
    if not adata.obs.index.is_integer():
        adata.obs["label_id"] = adata.obs.index
        adata.obs.index = pd.Categorical(adata.obs.index)
        adata.obs.index = adata.obs.index.codes
        adata.obs.index = adata.obs.index.astype(str)
        
    # turn obsm into a numpy array
    for k in adata.obsm_keys():
        adata.obsm[k] = np.array(adata.obsm[k])
    # print(adata.uns["spatial"][sample_id]["scalefactors"])

    spot_diameter_fullres = adata.uns["spatial"][sample_id]["scalefactors"][
        "spot_diameter_fullres"
    ]
    # hires_scalef = adata.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
    spot_coords = adata.obsm["spatial"]
    # print(spot_coords, X, Y, Z, C, T)
    assert adata.obs.shape[0] == spot_coords.shape[0]

    X = ome_md["X"]
    Y = ome_md["Y"]
    labelImg = np.zeros((int(Y), int(X)), dtype=np.uint16)
    # print(labelImg.shape)
    # print(np.max(spot_coords))
    for spId, (y, x) in zip(adata.obs.index, spot_coords):
        labelImg[disk((x, y), spot_diameter_fullres / 2)] = spId

    tf.imwrite(f"{stem}.tif", labelImg)


if __name__ == "__main__":
    fire.Fire(main)
