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


def main(stem, h5ad):
    adata = sc.read(h5ad)
    # sample_ids = np.unique(adata.obs[["sample"]])
    # print(sample_ids)
    # md = {}
    # for i in adata.uns["spatial"]:
        # sample_id = re.search(".*_count_\d+_(.+)_mm.*", i).group(1)
        # md[sample_id] = adata.uns["spatial"][i]["scalefactors"]
    # adata.obs["Y"] = adata.obsm["spatial"][:, 0]
    # adata.obs["X"] = adata.obsm["spatial"][:, 1]
    print(adata.uns["spatial"])
    print(adata.obsm["spatial"])


if __name__ == "__main__":
    fire.Fire(main)
