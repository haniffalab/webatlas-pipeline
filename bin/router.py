#!/usr/bin/env python3

import fire
import os

from process_h5ad import h5ad_to_zarr
from process_molecules import tsv_to_json
from process_spaceranger import spaceranger_to_zarr


def main(file_type, path, stem, args):
    out_file = None

    if file_type == "spaceranger":
        out_file = spaceranger_to_zarr(path=path, stem=stem, **args)
    elif file_type == "h5ad":
        out_file = h5ad_to_zarr(file=path, stem=stem, **args)
    elif file_type == "molecules":
        out_file = tsv_to_json(file=path, stem=stem, **args)

    return out_file


if __name__ == "__main__":
    fire.Fire(main)
