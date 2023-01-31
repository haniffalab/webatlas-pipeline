#!/usr/bin/env python3
"""
generate_label.py
====================================
Generates the label image from spatial data
"""

from __future__ import annotations
import os
import fire
import logging
import glob
import h5py
import zarr
import typing as T
import scanpy as sc
import numpy as np
from skimage.draw import disk, polygon
import tifffile as tf
import pandas as pd
from process_spaceranger import spaceranger_to_anndata


def visium(
    stem: str,
    file_path: str,
    shape: tuple[int, int] = None,
    obs_subset: tuple[int, T.Any] = None,
    sample_id: str = None,
) -> None:
    """This function writes a label image tif file with drawn labels according to an
    Anndata object with necessary metadata stored within `uns["spatial"]`.

    Args:
        stem (str): Prefix for the output image filename.
        file_path (str): Path to the h5ad file or spaceranger output directory.
        shape (tuple[int, int], optional): Output image shape. Defaults to None.
        obs_subset (tuple(str, T.Any), optional): Tuple containing an `obs` column name and one or more values
            to use to subset the AnnData object. Defaults to None.
        sample_id (str, optional): Sample ID string within the Anndata object. Defaults to None.
    """
    # sample_id = sample_id or Path(file_path).stem

    if os.path.isdir(file_path):
        adata = spaceranger_to_anndata(file_path)
        sample_id = sample_id or list(adata.uns["spatial"].keys())[0]
        if not shape:
            hires_shape = adata.uns["spatial"][sample_id]["images"]["hires"].shape
            scalef = adata.uns["spatial"][sample_id]["scalefactors"][
                "tissue_hires_scalef"
            ]
            shape = [int(hires_shape[0] / scalef), int(hires_shape[1] / scalef)]
    else:
        adata = sc.read(file_path)
        sample_id = sample_id or list(adata.uns["spatial"].keys())[0]

    # Subset adata by obs
    if obs_subset:
        obs_subset[1] = (
            [obs_subset[1]]
            if not isinstance(obs_subset[1], (list, tuple))
            else obs_subset[1]
        )
        adata = adata[adata.obs[obs_subset[0]].isin(obs_subset[1])]

    # check if index is numerical, if not reindex
    if not adata.obs.index.is_integer() and not (
        adata.obs.index.is_object() and all(adata.obs.index.str.isnumeric())
    ):
        adata.obs["label_id"] = adata.obs.index
        adata.obs.index = pd.Categorical(adata.obs.index)
        adata.obs.index = adata.obs.index.codes
        adata.obs.index = adata.obs.index.astype(str)

    # turn obsm into a numpy array
    for k in adata.obsm_keys():
        adata.obsm[k] = np.array(adata.obsm[k])

    spot_diameter_fullres = adata.uns["spatial"][sample_id]["scalefactors"][
        "spot_diameter_fullres"
    ]
    # hires_scalef = adata.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
    spot_coords = adata.obsm["spatial"]
    assert adata.obs.shape[0] == spot_coords.shape[0]

    label_img = np.zeros((shape[0], shape[1]), dtype=np.uint16)

    for spId, (y, x) in zip(adata.obs.index, spot_coords):
        label_img[disk((x, y), spot_diameter_fullres / 2)] = int(spId)

    tf.imwrite(f"{stem}-label.tif", label_img)

    return


def merscope(stem: str, path: str, shape: tuple[int, int], z_index: list[int] = [0]):
    """This function writes a label image tif file with drawn labels according
    to `cell_boundaries` data stored in MERSCOPE output directory

    Args:
        stem (str): Prefix for the output image filename.
        path (str): Path to the MERSCOPE output directory
        shape (tuple[int, int]): Output image shape.
        z_index (list[int], optional): Z indices to process. Defaults to [0].
    """

    z_index = [z_index] if not isinstance(z_index, (list, tuple)) else z_index

    tm = pd.read_csv(
        os.path.join(path, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep=" ",
        header=None,
        dtype=float,
    ).values

    fovs = [
        x
        for x in os.listdir(os.path.join(path, "cell_boundaries"))
        if x.endswith(".hdf5")
    ]

    for i, z in [(x, "zIndex_{}".format(x)) for x in z_index]:

        label_img = np.zeros((shape[0], shape[1]), dtype=np.uint32)

        for fov in fovs:
            with h5py.File(os.path.join(path, "cell_boundaries", fov)) as f:
                for cell_id in f["featuredata"].keys():
                    pol = f["featuredata"][cell_id][z]["p_0"]["coordinates"][0]

                    pol[:, 0] = pol[:, 0] * tm[0, 0] + tm[0, 2]
                    pol[:, 1] = pol[:, 1] * tm[1, 1] + tm[1, 2]

                    rr, cc = polygon(pol[:, 1], pol[:, 0])
                    label_img[rr - 1, cc - 1] = int(cell_id)

        logging.info(f"Writing label tif image {stem}-label_z{i} ...")
        tf.imwrite(
            f"{stem}-label_z{i}.tif" if len(z_index) > 1 else f"{stem}-label.tif",
            label_img,
        )

    return


def xenium(stem: str, path: str, shape: tuple[int, int], resolution: int = 0.2125):
    if os.path.isdir(path):
        cells_file = glob.glob(os.path.join(path, "*cells.zarr.zip"))[0]
    else:
        cells_file = path

    z = zarr.open(cells_file, "r")
    ids = z["cell_id"]
    pols = z["polygon_vertices"][1]

    label_img = np.zeros((shape[0], shape[1]), dtype=np.min_scalar_type(max(ids)))

    for id, pol in zip(ids, pols):
        pol = pol / resolution
        pol = np.array(list(map(list, pol.reshape(pol.shape[0] // 2, 2))))
        rr, cc = polygon(pol[:, 1], pol[:, 0])
        label_img[rr - 1, cc - 1] = int(id)

    tf.imwrite(f"{stem}-label.tif", label_img)

    return


def create_img(
    stem: str,
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

    if ref_img:
        tif_img = tf.TiffFile(ref_img)
        args["shape"] = tif_img.pages[0].shape[:2]

    if file_type == "visium":
        visium(stem, file_path, **args)
    elif file_type == "merscope":
        merscope(stem, file_path, **args)
    elif file_type == "xenium":
        xenium(stem, file_path, **args)


if __name__ == "__main__":
    fire.Fire(create_img)
