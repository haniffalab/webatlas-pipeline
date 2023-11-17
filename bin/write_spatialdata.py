#!/usr/bin/env python3
"""
write_spatialdata.py
====================================
Processes H5AD and images into SpatialData
"""

from __future__ import annotations
from typing import Union
import logging
import warnings
import fire
import tifffile as tf
import anndata as ad
import xarray as xr
import spatialdata as sd
from dask_image.imread import imread

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.INFO)


def read_image(path: str, is_label: bool = False):
    tif = tf.TiffFile(path)
    dims = list(tif.series[0].axes.lower()
                .replace("s", "c")
                .replace("i", "c")
               )
    image = imread(path).squeeze()
    imarray = xr.DataArray(image, dims=dims).chunk(chunks="auto")
    if is_label:
        return sd.models.Labels2DModel.parse(imarray)
    else:
        return sd.models.Image2DModel.parse(imarray)


def write_spatialdata(
    anndata_path: str,
    stem: str = "",
    raw_img_path: Union[str, list[str]] = [],
    label_img_path: Union[str, list[str]] = [],
) -> str:
    """This function takes an AnnData object and image files to
    write to a SpatialData objetc

    Args:
        path (str): Path to the h5ad file.
        stem (str, optional): Prefix for the output file. Defaults to "".
        raw_img_path (str, optional): Raw image to process. Defaults to None.
        label_img_path (str, optional): Label image to process. Defaults to None.

    Returns:
        str: Output SpatialData filename
    """
    if anndata_path.endswith(".h5ad"):
        adata = ad.read(anndata_path, backed=True)
    elif anndata_path.endswith(".zarr"):
        adata = ad.read_zarr(anndata_path)
    else:
        raise SystemError("Path to AnnData not .h5ad nor .zarr")

    OBS_IDX_NAME = "webatlas_index"
    adata.obs = adata.obs.reset_index(names=OBS_IDX_NAME).set_index(
        OBS_IDX_NAME, drop=False
    )  # have index as both index and column
    adata.obs[OBS_IDX_NAME] = adata.obs[OBS_IDX_NAME].astype(int)

    # ensure library_id in obs
    if "library_id" not in adata.obs:
        adata.obs["library_id"] = 0
        adata.obs["library_id"] = adata.obs["library_id"].astype("category")

    region_key = "library_id"
    region = adata.obs["library_id"].cat.categories.to_list()

    # `region_key`, `region` and `instance_key` are supposed to be optional
    # but enforced by code
    sd.models.TableModel.parse(
        adata, region_key=region_key, region=region, instance_key=OBS_IDX_NAME
    )

    sdata = sd.SpatialData(table=adata)

    if isinstance(raw_img_path, str):
        raw_img_path = [raw_img_path]
    for raw_img in raw_img_path:
        sdata.add_image("raw", read_image(raw_img))

    if isinstance(label_img_path, str):
        label_img_path = [label_img_path]
    for label_img in label_img_path:
        sdata.add_labels("label", read_image(label_img, is_label=True))

    zarr_file = f"{stem}-spatialdata.zarr"
    sdata.write(zarr_file)

    return zarr_file


if __name__ == "__main__":
    fire.Fire(write_spatialdata)
