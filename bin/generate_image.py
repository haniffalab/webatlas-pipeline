#!/usr/bin/env python3
"""
generate_image.py
====================================
Generates raw/label images from spatial data
"""

from __future__ import annotations
import fire
import typing as T
import tifffile as tf
from process_spaceranger import visium_label
from process_xenium import xenium_label
from process_merscope import merscope_label, merscope_raw


def create_img(
    stem: str,
    img_type: str,
    file_type: str,
    file_path: str,
    ref_img: str = None,
    args: dict[str, T.Any] = {},
) -> None:
    """This function calls the corresponding function
    to write a label image given the metadata provided.
    It also obtains the image shape of a reference image if specified.

    Args:
        stem (str): Prefix for the output image filename.
        file_type (str): Type of file containing the metadata from which to
            generate the label image.
        file_path (str): Path to the metadata file.
        ref_img (str, optional): Path to reference image from which to get the
            shape for the label image. Defaults to None.
        args (dict[str,T.Any], optional): Args to be passed to the appropriate processing function.
            Defaults to {}.
    """
    print(args)
    if ref_img:
        tif_img = tf.TiffFile(ref_img)
        args["shape"] = tif_img.pages[0].shape[:2]

    if img_type == "label":
        if file_type == "visium":
            visium_label(stem, file_path, **args)
        elif file_type == "merscope":
            merscope_label(stem, file_path, **args)
        elif file_type == "xenium":
            xenium_label(stem, file_path, **args)
    elif img_type == "raw":
        if file_type == "merscope":
            merscope_raw(stem, file_path, **args)


if __name__ == "__main__":
    fire.Fire(create_img)
