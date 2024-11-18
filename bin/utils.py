from __future__ import annotations
import tifffile as tf


def tiff_image_size(
    tiff_file: str,
):
    tif_img = tf.TiffFile(tiff_file)
    return tif_img.pages[0].shape[:2]


def visium_image_size(
    adata: sc.AnnData,
):
    sample = list(adata.uns["spatial"].keys())[0]
    hiresfactor = adata.uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"]
    n, m, _ = adata.uns["spatial"][sample]["images"]["hires"].shape
    m = int(m / hiresfactor)
    n = int(n / hiresfactor)
    return [m, n]
