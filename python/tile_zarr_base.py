#!/usr/bin/env python3

import numpy as np
import dask.array as da
from numcodecs import Zlib

from pathlib import Path
import argparse


def tile_zarr(
    zarr_pyramid_base,
    max_level=None,
    compressor=None,
    dtype=None,
):
    img = da.from_zarr(zarr_pyramid_base)
    pyramid_path = Path(zarr_pyramid_base).parent
    chunks = img.chunksize
    tile_size = int(img.chunksize[-1])

    if dtype is None:
        dtype = img.dtype
    if max_level is None:
        # create all levels up to 512 x 512
        max_level = int(np.ceil(np.log2(np.maximum(img.shape[1], img.shape[2])))) - 9
    if compressor is None:
        compressor = Zlib(level=1)

    for i in range(1, max_level):
        img = da.coarsen(np.mean, img, {1: 2, 2: 2}, trim_excess=True)

        # Edge Case: Need to pad smallest thumbnail sometimes.
        #
        # If a dimension of an array is larger than TILE_SIZE,
        # zarr will respect the chunk size requested an automatically
        # pad with zeros in the store. However, if an array dimension
        # is smaller than the tile size, `da.to_zarr` will change the
        # chunking and not pad with zeros. We need sometimes need to pad
        # for the smallest tiles because x and y might not be square.
        if img.shape[1] < tile_size:
            img = da.pad(
                img,
                ((0, 0), (0, tile_size - img.shape[1]), (0, 0)),
                "constant",
            )

        if img.shape[2] < tile_size:
            img = da.pad(
                img,
                ((0, 0), (0, 0), (0, tile_size - img.shape[2])),
                "constant",
            )

        # Define pyramid level path
        out_path = str(pyramid_path / str(i))

        # Write to zarr store
        img.astype(dtype).rechunk(chunks).to_zarr(out_path, compressor=compressor)

        # Read from last store so dask doesn't need to re-compute
        # task graph starting at base.
        img = da.from_zarr(out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create zarr tile pyramid from base resolution zarr"
    )
    parser.add_argument(
        "--zarr_pyramid_base", required=True, help="zarr store with base image"
    )
    args = parser.parse_args()

    tile_zarr(args.zarr_pyramid_base)
