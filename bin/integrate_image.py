#!/usr/bin/env python3

import fire
import os
import tifffile as tf
import numpy as np
import logging
from ome_zarr.reader import Reader
from ome_zarr.io import parse_url
import shutil
from ome_zarr.writer import write_multiscale, write_image
import zarr


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
    label = nodes[0].data[0]
    reindexed_label = add_offset(label.squeeze(), offset)
    os.makedirs(f"{out_filename}/OME", exist_ok=True)
    store = parse_url("./", mode="w").store
    temp_group = zarr.group(store=store).create_group(out_filename)
    # write_multiscale(
    #     np.array(reindexed_label), temp_group, axes="yx"
    # )  # this should be used, but need to bump ome-zarr version
    write_image(np.array(reindexed_label), temp_group, axes="yx")
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
