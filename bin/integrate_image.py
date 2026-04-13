#!/usr/bin/env python3

import logging
import os
import shutil

import fire
import numpy as np
import tifffile as tf
import zarr
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
from ome_zarr.writer import write_multiscale


def add_offset(label, offset: int):
    max_id = np.max(label)
    mask = label != 0
    reindexed_label = (label + offset) * mask
    reindexed_max_id = np.max(reindexed_label)
    logging.info(max_id, reindexed_max_id)
    return reindexed_label


def reindex_label(label_image: str, offset: int, out_filename: str) -> None:
    label = tf.imread(label_image).astype(np.int32)
    reindexed_label = add_offset(label, offset)
    tf.imwrite(out_filename, reindexed_label)


def reindex_label_zarr(label_image_path: str, offset: int, out_filename: str) -> None:
    binary_path = (
        label_image_path
        if label_image_path.endswith("/0")
        else os.path.join(label_image_path, "0")
    )
    reader = Reader(parse_url(binary_path))
    nodes = list(reader())
    labels = nodes[0].data
    reindexed_labels = [add_offset(x, offset) for x in labels]

    zarr_format = (
        3 if os.path.exists(os.path.join(label_image_path, "zarr.json")) else 2
    )

    store = parse_url(out_filename, mode="w").store
    tmp_group = zarr.open_group(store=store, mode="w", zarr_format=zarr_format)
    write_multiscale(
        reindexed_labels,
        tmp_group,
        compute=True,
        storage_options=dict(dimension_separator="/"),
    )
    zarr.consolidate_metadata(out_filename)
    os.makedirs(f"{out_filename}/OME", exist_ok=True)
    shutil.copy(
        label_image_path + "/OME/METADATA.ome.xml",
        f"{out_filename}/OME/METADATA.ome.xml",
    )


def process_image(label_image_path: str, **kwargs) -> None:
    ext = os.path.splitext(label_image_path)[-1]
    if ext.lower() in [".tif", ".tiff"]:
        reindex_label(label_image_path, **kwargs)
    elif ext.lower() in [".zarr"]:
        reindex_label_zarr(label_image_path, **kwargs)
    return


if __name__ == "__main__":
    fire.Fire(process_image)
