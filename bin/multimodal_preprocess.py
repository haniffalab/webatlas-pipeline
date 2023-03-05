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


def reindex_label(label_image:str, start_at: int, out_name:str) -> None:
    lab = tf.imread(label_image).astype(np.int32)
    max_id = np.max(lab)
    mask = lab!=0
    reindexed_lab = (lab + start_at) * mask
    reindexed_max_id = np.max(lab)
    logging.info(max_id, reindexed_max_id)
    tf.imwrite(out_name, reindexed_lab)

if __name__ == "__main__":
    fire.Fire(reindex_label)
