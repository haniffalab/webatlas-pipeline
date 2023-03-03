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
from process_xenium import xenium_to_zarr


def main(file_type: str, path: str, stem: str, args: dict[str, T.Any] = {}) -> str:
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

    if file_type == "spaceranger":
        out_file = spaceranger_to_zarr(path=path, stem=stem, **args)
    elif file_type == "h5ad":
        out_file = h5ad_to_zarr(file=path, stem=stem, **args)
    elif file_type == "molecules":
        out_file = tsv_to_json(file=path, stem=stem, **args)
    elif file_type == "xenium":
        out_file = xenium_to_zarr(file=path, stem=stem, **args)

    return out_file


if __name__ == "__main__":
    fire.Fire(main)
