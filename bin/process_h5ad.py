#!/usr/bin/env python3
"""
process_h5ad.py
====================================
Processes H5AD files into AnnData-Zarr
"""

from __future__ import annotations

import logging
import os
import typing as T
import warnings

import anndata as ad
import fire
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import zarr
from constants.suffixes import ANNDATA_ZARR_SUFFIX
from scipy.sparse import csc_matrix, csr_matrix, spmatrix
from tqdm import tqdm

from utils import visium_image_size

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.INFO)
ad.settings.zarr_write_format = 3
ad.settings.allow_write_nullable_strings = False
zarr.config.set({"async.concurrency": 1})


def h5ad_to_zarr(
    path: str = None,
    stem: str = "",
    adata: ad.AnnData = None,
    chunk_size: int = 10,
    batch_processing: bool = False,
    batch_size: int = 10000,
    consolidate_metadata: bool = True,
    convert_strings_to_categoricals: bool = False,
    chunks_per_shard: int = 100,
    append: bool = False,
    **kwargs,
) -> str:
    """This function takes an AnnData object or path to an h5ad file,
    ensures data is of an appropriate data type for Vitessce
    and writes the object to Zarr.

    Args:
        path (str, optional): Path to the h5ad file. Defaults to None.
        stem (str, optional): Prefix for the output file. Defaults to "".
        adata (AnnData, optional): AnnData object to process. Supersedes `path`.
            Defaults to None.
        chunk_size (int, optional): Output Zarr column chunk size. Defaults to 10.
        batch_processing (bool, optional): If the expression matrix will be written
            to Zarr incrementally. Use to avoid loading the whole AnnData into memory.
            Defaults to False.
        batch_size (int, optional): The amount of rows (if matrix is in CSR format)
            or columns (if matrix is dense or in CSC format) of the expression matrix
            to process at a time when batch processing. Defaults to 10000.
        consolidate_metadata (bool, optional): Whether to consolidate Zarr metadata
            after writing. Defaults to True.
        convert_strings_to_categoricals (bool, optional): Whether to convert string
            columns to categorical when writing to Zarr. Defaults to False.
        chunks_per_shard (int, optional): Number of chunks per shard for Zarr
            sharded arrays. Defaults to 100.
        append (bool, optional): If True, append to existing Zarr array instead of
            creating new. Defaults to False.

    Raises:
        SystemError: If `batch_processing` is True and the matrix contains
            an `indptr` key but the matrix is not in scipy.sparse.csr_matrix
            nor scipy.sparse.csc_matrix format

    Returns:
        str: Output Zarr directory path
    """

    zarr_dir = (
        f"{stem}-{ANNDATA_ZARR_SUFFIX}"
        if not stem.endswith("-" + os.path.splitext(ANNDATA_ZARR_SUFFIX)[0])
        else f"{stem}{os.path.splitext(ANNDATA_ZARR_SUFFIX)[1]}"
    )

    if adata is not None:
        adata = preprocess_anndata(adata, **kwargs)

        m = len(adata.obs)

        X = adata.X
        del adata.X

        adata.write_zarr(
            zarr_dir,
            convert_strings_to_categoricals=convert_strings_to_categoricals,
        )

        # @TODO: support batch processing to avoid .toarray() of whole matrix
        # matrix sparse to dense
        if isinstance(X, spmatrix):
            # use toarray() as it generates a ndarray
            # instead of todense() which generates a matrix
            X = X.toarray()

        if chunks_per_shard:
            shards = (m, chunk_size * chunks_per_shard)
        else:
            shards = None

        zarr.create_array(
            zarr_dir,
            name="X",
            data=X,
            chunks=(m, chunk_size),
            shards=shards,
        )

        if consolidate_metadata:
            zarr.consolidate_metadata(zarr_dir)

        logging.info("Wrote AnnData object to {}".format(zarr_dir))

        return zarr_dir

    with h5py.File(path, "r") as f:
        adata = ad.AnnData(
            obs=ad._io.h5ad.read_elem(f["obs"]) if "obs" in f else None,
            var=ad._io.h5ad.read_elem(f["var"]) if "var" in f else None,
            obsm=ad._io.h5ad.read_elem(f["obsm"]) if "obsm" in f else None,
            obsp=ad._io.h5ad.read_elem(f["obsp"]) if "obsp" in f else None,
            varm=ad._io.h5ad.read_elem(f["varm"]) if "varm" in f else None,
            varp=ad._io.h5ad.read_elem(f["varp"]) if "varp" in f else None,
            uns=(ad._io.h5ad.read_elem(f["uns"]) if "uns" in f else None),
        )

        m = len(adata.obs)
        n = len(adata.var)

        logging.info("Matrix shape: {} x {}".format(m, n))

        is_sparse = isinstance(f["X"], h5py.Group) and "indptr" in f["X"].keys()
        if is_sparse:
            if len(f["X"]["indptr"]) - 1 == m:
                is_csc = False
            elif len(f["X"]["indptr"]) - 1 == n:
                is_csc = True
            else:
                raise SystemError("Error identifying sparse matrix format")

    adata = preprocess_anndata(adata, **kwargs)

    # Write anndata tables to zarr
    logging.info("Writing anndata tables")
    adata.write_zarr(
        zarr_dir, convert_strings_to_categoricals=convert_strings_to_categoricals
    )

    if os.path.exists(os.path.join(zarr_dir, ".zmetadata")):
        os.remove(os.path.join(zarr_dir, ".zmetadata"))

    if batch_processing:
        batch_size = max(1, batch_size)
        logging.info("Batch processing with batch size {}".format(batch_size))
    else:
        batch_size = m if is_sparse and not is_csc else n

    if chunks_per_shard:
        shards = (m, chunk_size * chunks_per_shard)
    else:
        shards = None

    if append:
        z = zarr.create_array(
            zarr_dir,
            name="X",
            shape=(m, 0) if not is_sparse else (m if is_csc else 0, 0 if is_csc else n),
            chunks=(m, chunk_size),
            dtype="float32",
            shards=shards,
        )
    else:
        z = zarr.create_array(
            zarr_dir,
            name="X",
            shape=(m, n),
            chunks=(m, chunk_size),
            dtype="float32",
            shards=shards,
            overwrite=True,
        )

    logging.info("Chunk shape: {}. Shards: {}".format(z.chunks, z.shards))

    if is_sparse:
        batch_process_sparse(
            path,
            zarr_dir=zarr_dir,
            m=m,
            n=n,
            batch_size=batch_size,
            chunk_size=chunk_size,
            shards=shards,
            is_csc=is_csc,
            append=append,
        )
    else:
        batch_process_array(
            path,
            zarr_dir=zarr_dir,
            m=m,
            n=n,
            batch_size=batch_size,
            chunk_size=chunk_size,
            shards=shards,
            append=append,
        )

    if consolidate_metadata:
        zarr.consolidate_metadata(zarr_dir)

    logging.info("Wrote AnnData object to {}".format(zarr_dir))

    return zarr_dir


def reindex_anndata_obs(adata: ad.AnnData) -> ad.AnnData:
    # check if index is numerical, if not reindex
    if not adata.obs.index.is_integer() and not (
        adata.obs.index.is_object() and all(adata.obs.index.str.isnumeric())
    ):
        IDX_NAME = "label_id"
        if IDX_NAME in adata.obs:
            adata.obs.rename(columns={IDX_NAME: f"_{IDX_NAME}"})
        adata.obs = adata.obs.reset_index(names=IDX_NAME)
        adata.obs.index = (
            pd.Categorical(adata.obs[IDX_NAME]).codes + 1
        )  # avoid 0's for possible label images
        adata.uns["webatlas_reindexed"] = True
    adata.obs.index = adata.obs.index.astype(str)

    return adata


def subset_anndata(
    adata: ad.AnnData,
    obs_subset: tuple[str, T.Any] = None,
    var_subset: tuple[str, T.Any] = None,
    sample: str = None,
) -> ad.AnnData:
    # Subset adata by obs
    if obs_subset:
        logging.info(f"Subsetting AnnData by {obs_subset[0]}")
        obs_subset[1] = (
            [obs_subset[1]]
            if not isinstance(obs_subset[1], (list, tuple))
            else obs_subset[1]
        )
        adata = adata[adata.obs[obs_subset[0]].isin(obs_subset[1])]

    # Subset adata by var
    if var_subset:
        logging.info(f"Subsetting AnnData by {var_subset[0]}")
        var_subset[1] = (
            [var_subset[1]]
            if not isinstance(var_subset[1], (list, tuple))
            else var_subset[1]
        )
        adata = adata[:, adata.var[var_subset[0]].isin(var_subset[1])]

    # Remove other samples' spatial data
    if sample and "spatial" in adata.uns:
        spatial_samples = list(adata.uns["spatial"].keys())
        for spatial_sample in spatial_samples:
            if spatial_sample != sample:
                del adata.uns["spatial"][spatial_sample]

    return adata


def rotate_anndata(
    adata: ad.AnnData,
    shape: tuple[int, int],
    degrees: T.Literal[90, 180, 270],
) -> ad.AnnData:
    """
    Counterclockwise rotate the spatial coordinates and images in an AnnData object
    """
    if degrees not in [90, 180, 270]:
        raise SystemError("Invalid rotation degrees. Must be 90, 180, or 270.")

    logging.info(f"Rotating spatial coordinates and images by {degrees} degrees")

    m, n = shape[0], shape[1]

    for spatial_key in ["spatial", "X_spatial"]:
        if spatial_key in adata.obsm:
            rot_spatial = []
            for [x, y] in adata.obsm[spatial_key]:
                if degrees == 90:
                    rot_spatial.append([y, m - x])
                elif degrees == 180:
                    rot_spatial.append([m - x, n - y])
                elif degrees == 270:
                    rot_spatial.append([n - y, x])

            adata.obsm[spatial_key] = np.array(rot_spatial)

    adata.uns["webatlas_rotation"] = degrees

    return adata


def rescale_spatial(adata: ad.AnnData, factor: float) -> ad.AnnData:
    """
    Rescale spatial coordinates in an AnnData object
    """

    logging.info(f"Rescaling spatial coordinates by {factor}")

    for spatial_key in ["spatial", "X_spatial"]:
        if spatial_key in adata.obsm:
            adata.obsm[spatial_key] = adata.obsm[spatial_key] * factor

    adata.uns["webatlas_rescale"] = factor

    return adata


def preprocess_anndata(
    adata: ad.AnnData,
    compute_embeddings: bool = False,
    var_index: str = None,
    obs_subset: tuple[str, T.Any] = None,
    var_subset: tuple[str, T.Any] = None,
    spatial_shape: tuple[int, int] = None,
    rotate_degrees: T.Literal[90, 180, 270] = None,
    rescale_factor: float = None,
    sample: str = None,
    **kwargs,
):
    """This function preprocesses an AnnData object, ensuring correct dtypes
        for zarr conversion

    Args:
        adata (AnnData): AnnData object to preprocess.
        compute_embeddings (bool, optional): If `X_umap` and `X_pca` embeddings
            will be computed. Defaults to False.
        var_index (str, optional): Alternative `var` column name with `var` names
            to be used in the visualization. Defaults to None.
        obs_subset (tuple(str, T.Any), optional): Tuple containing an `obs` column name
            and one or more values to use to subset the AnnData object.
            Defaults to None.
        var_subset (tuple(str, T.Any), optional): Tuple containing a `var` column name
            and one or more values to use to subset the AnnData object.
            Defaults to None.
        spatial_shape (tuple[int, int], optional): Shape (height, width)
            for spatial rotation. Defaults to None.
        rotate_degrees (Literal[90, 180, 270], optional): Degrees to rotate
            spatial coordinates. Defaults to None.
        rescale_factor (float, optional): Factor to rescale spatial coordinates.
            Defaults to None.
        sample (str, optional): Sample name to filter spatial data. Defaults to None.

    Returns:
        AnnData: Preprocessed AnnData object ready for Zarr conversion.
    """

    adata = subset_anndata(
        adata, obs_subset=obs_subset, var_subset=var_subset, sample=sample
    )

    if rescale_factor:
        adata = rescale_spatial(adata, rescale_factor)

    if rotate_degrees:
        if not spatial_shape:
            try:
                spatial_shape = visium_image_size(adata)
            except Exception:
                raise SystemError("Must provide spatial shape to rotate spatial data.")
        adata = rotate_anndata(adata, spatial_shape, rotate_degrees)

    # ensure data types for var
    for col_name, col_data in adata.var.reset_index().items():
        mod = False
        if col_data.dtype == "string":
            col_data = col_data.astype("str")
            mod = True
        elif col_data.dtype == "category" and col_data.cat.categories.dtype == "string":
            col_data = col_data.astype("str").astype("category")
            mod = True

        if mod:
            if col_name == adata.var.index.name:
                adata.var.index = col_data
            else:
                adata.var[col_name] = col_data

    # reindex var with a specified column
    if var_index and var_index in adata.var:
        try:
            adata.var.reset_index(inplace=True)
        except ValueError:
            logging.warning(
                "Column already exists when trying to reset var index. Dropping index."
            )
            adata.var.reset_index(inplace=True, drop=True)
        adata.var.set_index(var_index, inplace=True)
        adata.var.index = adata.var.index.astype(str)
    adata.var_names_make_unique()

    # turn obsm into a numpy array
    for k in adata.obsm_keys():
        adata.obsm[k] = np.array(adata.obsm[k])

    adata = reindex_anndata_obs(adata)

    # compute embeddings if not already stored in object
    if compute_embeddings:
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata)
        if "X_umap" not in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

    # ensure data types for obs
    for col_name, col_data in adata.obs.reset_index().items():
        mod = False
        if col_data.dtype == "string":
            col_data = col_data.astype("str")
            mod = True
        elif col_data.dtype == "category" and col_data.cat.categories.dtype == "string":
            col_data = col_data.astype("str").astype("category")
            mod = True
        elif col_data.dtype in ["int8", "int64"]:
            col_data = col_data.astype("int32")
            mod = True
        elif col_data.dtype == "bool":
            col_data = col_data.astype(str).astype("category")
            mod = True

        if mod:
            if col_name == adata.obs.index.name:
                adata.obs.index = col_data
            else:
                adata.obs[col_name] = col_data

    # ensure data types for obsm
    for col in adata.obsm:
        if type(adata.obsm[col]).__name__ in ["DataFrame", "Series"]:
            adata.obsm[col] = adata.obsm[col].to_numpy()
        if adata.obsm[col].dtype in ["int8", "int64"] or col == "spatial":
            adata.obsm[col] = adata.obsm[col].astype("int32")

    # remove unnecessary data
    del adata.raw

    return adata


def batch_process_sparse(
    file: str,
    zarr_dir: str,
    m: int,
    n: int,
    batch_size: int,
    chunk_size: int,
    shards: str = None,
    is_csc: bool = False,
    append: bool = False,
) -> None:
    """Function to incrementally load and write a sparse matrix to Zarr

    Args:
        file (str): Path to h5ad file.
        zarr_dir (str): Path to output Zarr directory.
        m (int): Number of rows in the matrix.
        n (int): Number of columns in the matrix.
        batch_size (int): Number of rows/columns to load and write at a time.
        chunk_size (int): Output Zarr column chunk size.
        shards (str, optional): Zarr shard configuration. Defaults to None.
        is_csc (bool, optional): If matrix is in CSC format instead of CSR format.
            Defaults to False.
        append (bool, optional): If True, append to existing array. Defaults to False.
    """

    logging.info("Processing sparse {} matrix".format("CSC" if is_csc else "CSR"))

    z = zarr.open_array(zarr_dir, path="X", mode="a")

    compressed_dim = n if is_csc else m

    with h5py.File(file, "r") as f:
        indptr = f["X"]["indptr"][:]
        batch_size = compressed_dim if batch_size > compressed_dim else batch_size
        for i in tqdm(range((compressed_dim + batch_size - 1) // batch_size)):
            j = i * batch_size
            if j + batch_size > compressed_dim:
                batch_size = compressed_dim - j
            k = j + batch_size

            indices = f["X"]["indices"][indptr[j] : indptr[k]]
            data = f["X"]["data"][indptr[j] : indptr[k]]

            if is_csc:
                matrix = csc_matrix(
                    (data, indices, indptr[j : k + 1] - indptr[j]),
                    shape=(m, 1 * batch_size),
                )
            else:
                matrix = csr_matrix(
                    (data, indices, indptr[j : k + 1] - indptr[j]),
                    shape=(1 * batch_size, n),
                )

            if append:
                z.append(matrix.toarray(), axis=(1 * is_csc))
            else:
                if is_csc:
                    z[:, j:k] = matrix.toarray()
                else:
                    z[j:k, :] = matrix.toarray()
    return


def batch_process_array(
    file: str,
    zarr_dir: str,
    m: int,
    n: int,
    batch_size: int,
    chunk_size: int,
    shards: str = None,
    append: bool = False,
) -> None:
    """Function to incrementally load and write a dense matrix to Zarr

    Args:
        file (str): Path to h5ad file.
        zarr_dir (str): Path to output Zarr directory.
        m (int): Number of rows in the matrix.
        n (int): Number of columns in the matrix.
        batch_size (int): Number of columns to load and write at a time.
        chunk_size (int): Output Zarr column chunk size.
        shards (str, optional): Zarr shard configuration. Defaults to None.
        append (bool, optional): If True, append to existing array. Defaults to False.
    """

    logging.info("Processing dense matrix")

    z = zarr.open_array(zarr_dir, path="X", mode="a")

    with h5py.File(file, "r") as f:
        for i in tqdm(range((n + batch_size - 1) // batch_size)):
            j = i * batch_size
            if j + batch_size > n:
                batch_size = n - j
            k = j + batch_size

            matrix = f["X"][:, j:k]

            if append:
                z.append(matrix, axis=1)
            else:
                z[:, j:k] = matrix
    return


if __name__ == "__main__":
    fire.Fire(h5ad_to_zarr)
