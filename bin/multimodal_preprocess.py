#!/usr/bin/env python3

"""
multimodal_preprocess.py
====================================
Preprocesses integrated multimodal datasets to ensure optimal browsing performance
"""
import fire
import tifffile as tf
import numpy as np
import logging


def reindex_label(label_image: str, offset: int, out_name: str) -> None:
    label = tf.imread(label_image).astype(np.int32)
    max_id = np.max(label)
    mask = label != 0
    reindexed_label = (label + offset) * mask
    reindexed_max_id = np.max(label)
    logging.info(max_id, reindexed_max_id)
    tf.imwrite(out_name, reindexed_label)


if __name__ == "__main__":
    fire.Fire(reindex_label)
