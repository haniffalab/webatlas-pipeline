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
import pickle
import json
from h5ad_2_json import cell_sets_json, get_clusters


def main(dict_file, cells_file, cell_sets_file, matrix_file):
    with open(dict_file, 'rb') as handle:
        md = pickle.load(handle)

    json.dump(md, open(cells_file, "w"), indent=1)
    json.dump(cell_sets_json(md), open(cell_sets_file, "w"), indent=1)
    json.dump(get_clusters(md), open(matrix_file, "w"), indent=1)


if __name__ == "__main__":
    fire.Fire(main)
