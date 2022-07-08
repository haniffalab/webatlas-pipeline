#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import logging
import fire
from zarr import consolidate_metadata
from pathlib import Path
import zarr


def main(file_in):
    stem = Path(file_in).stem
    logging.info(stem)
    consolidated = consolidate_metadata(file_in)
    logging.info(consolidated.info)


if __name__ == "__main__":
    fire.Fire(main)
