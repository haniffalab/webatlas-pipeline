#!/usr/bin/env python3

import fire
import os
import tifffile as tf
import numpy as np
import logging


def reindex_label(label_image: str, offset: int, out_filename: str) -> None:
    label = tf.imread(label_image).astype(np.int32)
    max_id = np.max(label)
    mask = label != 0
    reindexed_label = (label + offset) * mask
    reindexed_max_id = np.max(label)
    logging.info(max_id, reindexed_max_id)
    tf.imwrite(out_filename, reindexed_label)


def process_image(label_image: str, **kwargs) -> None:

    ext = os.path.splitext(label_image)[-1]
    if ext.lower() in [".tif", ".tiff"]:
        reindex_label(**kwargs)
    # TODO: add option for zarr

    return


if __name__ == "__main__":
    fire.Fire(process_image)
