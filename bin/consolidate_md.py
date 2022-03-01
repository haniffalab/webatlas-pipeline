#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
from zarr import consolidate_metadata, save
from pathlib import Path
import zarr


def main(file_in):
    stem = Path(file_in).stem
    print(stem)
    # store = zarr.open(file_in)
    # print(store)
    consolidated = consolidate_metadata(file_in)
    print(consolidated.info)
    # save(f'./test_{stem}.zarr', consolidated)


if __name__ == "__main__":
    fire.Fire(main)
