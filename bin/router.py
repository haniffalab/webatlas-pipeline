#!/usr/bin/env python3
"""
router.py
====================================
Calls file processing functions
"""

from __future__ import annotations
import fire
import typing as T

from process_h5ad import h5ad_to_zarr
from process_molecules import tsv_to_json
from process_spaceranger import spaceranger_to_zarr
from process_merscope import merscope_to_zarr
from process_xenium import xenium_to_zarr


def process(file_type: str, path: str, stem: str, args: dict[str, T.Any] = {}) -> str:
    """Function that calls the appropriate processing function
    for the input file according to its type

    Args:
        file_type (str): Type of file to process
        path (str): Path to file to process
        stem (str): Prefix for output files
        args (dict[str,T.Any], optional): Args to be passed to the appropriate processing function.
            Defaults to {}.

    Returns:
        str: Output filename
    """
    out_file = None

    func_dict = {
        "spaceranger": spaceranger_to_zarr,
        "xenium": xenium_to_zarr,
        "merscope": merscope_to_zarr,
        "h5ad": h5ad_to_zarr,
        "molecules": tsv_to_json
    }

    out_file = func_dict[file_type](path=path, stem=stem, **args)

    return out_file


if __name__ == "__main__":
    fire.Fire(process)
