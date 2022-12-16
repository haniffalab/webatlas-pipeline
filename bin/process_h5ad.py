#!/usr/bin/env python3

import fire
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import h5py
import zarr
import logging
import warnings
from scipy.sparse import spmatrix, csr_matrix, csc_matrix

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.INFO)

SUFFIX = "anndata.zarr"


def h5ad_to_zarr(
    file: str = None,
    stem: str = "",
    adata: ad.AnnData = None,
    compute_embeddings: bool = False,
    chunk_size: int = 10,
    var_index: str = None,
    batch_processing: bool = False,
    batch_size: int = 10000,
) -> str:
    """This function takes an AnnData object or path to an h5ad file,
    ensures data is of an appropriate data type for Vitessce
    and writes the object to Zarr.

    Args:
        file (str, optional): Path to the h5ad file. Defaults to None.
        stem (str, optional): Prefix for the output file. Defaults to "".
        adata (AnnData, optional): AnnData object to process. Supersedes `file`.
            Defaults to None.
        compute_embeddings (bool, optional): If `X_umap` and `X_pca` embeddings will be computed.
            Defaults to False.
        chunk_size (int, optional): Output Zarr column chunk size. Defaults to 10.
        var_index (str, optional): Alternative `var` column name with `var` names
            to be used in the visualization. Defaults to None.
        batch_processing (bool, optional): If the expression matrix will be written to Zarr incrementally.
            Use to avoid loading the whole AnnData into memory. Defaults to False.
        batch_size (int, optional): The amount of rows (if matrix is in CSR format)
            or columns (if matrix is dense or in CSC format) of the expression matrix
            to process at a time when batch processing. Defaults to 10000.

    Raises:
        SystemError: If `batch_processing` is True and the matrix contains an `indptr` key
        but the matrix is not in scipy.sparse.csr_matrix nor scipy.sparse.csc_matrix format

    Returns:
        str: Output Zarr filename
    """

    if not adata:
        if not batch_processing:
            adata = sc.read(file)
        else:
            batch_size = max(1, batch_size)
            logging.info("Batch processing with batch size {}".format(batch_size))
            with h5py.File(file, "r") as f:
                # Load everything but X and layers
                adata = ad.AnnData(
                    obs=ad._io.h5ad.read_elem(f["obs"]) if "obs" in f else None,
                    var=ad._io.h5ad.read_elem(f["var"]) if "var" in f else None,
                    obsm=ad._io.h5ad.read_elem(f["obsm"]) if "obsm" in f else None,
                    obsp=ad._io.h5ad.read_elem(f["obsp"]) if "obsp" in f else None,
                    varm=ad._io.h5ad.read_elem(f["varm"]) if "varm" in f else None,
                    varp=ad._io.h5ad.read_elem(f["varp"]) if "varp" in f else None,
                    uns=ad._io.h5ad.read_elem(f["uns"]) if "uns" in f else None,
                )

    # reindex var with a specified column
    if var_index and var_index in adata.var:
        adata.var.reset_index(inplace=True)
        adata.var.set_index(var_index, inplace=True)
        adata.var.index = adata.var.index.astype(str)
    adata.var_names_make_unique()

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

    # compute embeddings if not already stored in object
    if compute_embeddings:
        if not "X_pca" in adata.obsm:
            sc.tl.pca(adata)
        if not "X_umap" in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

    for col in adata.obs:
        # anndata >= 0.8.0
        # if data type is categorical vitessce will throw "path obs/X contains a group" and won"t find .zarray
        # if adata.obs[col].dtype == "category":
        #     adata.obs[col] = adata.obs[col].cat.codes
        if adata.obs[col].dtype in ["int8", "int64"]:
            adata.obs[col] = adata.obs[col].astype("int32")
        if adata.obs[col].dtype == "bool":
            adata.obs[col] = adata.obs[col].astype(str).astype("category")

    for col in adata.obsm:
        if type(adata.obsm[col]).__name__ in ["DataFrame", "Series"]:
            adata.obsm[col] = adata.obsm[col].to_numpy()
        if adata.obsm[col].dtype in ["int8", "int64"] or col == "spatial":
            adata.obsm[col] = adata.obsm[col].astype("int32")

    # remove unnecessary data
    del adata.raw

    zarr_file = f"{stem}_{SUFFIX}"

    if not batch_processing:
        # matrix sparse to dense
        if isinstance(adata.X, spmatrix):
            # use toarray() as it generates a ndarray, instead of todense() which generates a matrix
            adata.X = adata.X.toarray()

        adata.write_zarr(zarr_file, [adata.shape[0], chunk_size])
    elif file and batch_processing:
        adata.write_zarr(zarr_file)

        m = len(adata.obs)
        n = len(adata.var)

        with h5py.File(file, "r") as f:
            if isinstance(f["X"], h5py.Group) and "indptr" in f["X"].keys():
                if len(f["X"]["indptr"]) - 1 == m:
                    logging.info("Batch processing sparse CSR matrix...")
                    batch_process_sparse(file, zarr_file, m, n, batch_size, chunk_size)
                elif len(f["X"]["indptr"]) - 1 == n:
                    logging.info("Batch processing sparse CSC matrix...")
                    batch_process_sparse(
                        file, zarr_file, m, n, batch_size, chunk_size, is_csc=True
                    )
                else:
                    raise SystemError("Error identifying sparse matrix format")
            else:
                logging.info("Batch processing dense matrix...")
                batch_process_array(file, zarr_file, m, n, batch_size, chunk_size)

    return zarr_file


def batch_process_sparse(
    file: str,
    zarr_file: str,
    m: int,
    n: int,
    batch_size: int,
    chunk_size: int,
    is_csc: bool = False,
) -> None:
    """Function to incrementally load and write a sparse matrix to Zarr

    Args:
        file (str): Path to h5ad file
        zarr_file (str): Path to output Zarr file
        m (int): Number of rows in the matrix
        n (int): Number of columns in the matrix
        batch_size (int): Number of rows/columns to load and write at a time
        chunk_size (int): Output Zarr column chunk size
        is_csc (bool, optional): If matrix is in CSC format instead of CSR format.
            Defaults to False.
    """
    l = n if is_csc else m
    z = zarr.open(zarr_file, mode="a")
    z["X"] = zarr.empty(
        (m if is_csc else 0, 0 if is_csc else n),
        chunks=(m, chunk_size),
        dtype="float32",
    )
    with h5py.File(file, "r") as f:
        indptr = f["X"]["indptr"][:]
        batch_size = l if batch_size > l else batch_size
        for i in range(l // batch_size + 1):
            j = i * batch_size
            if j + batch_size > l:
                batch_size = l - j
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

            z["X"].append(matrix.toarray(), axis=(1 * is_csc))
    return


def batch_process_array(
    file: str, zarr_file: str, m: int, n: int, batch_size: int, chunk_size: int
) -> None:
    """Function to incrementally load and write a dense matrix to Zarr

    Args:
        file (str): Path to h5ad file
        zarr_file (str): Path to output Zarr file
        m (int): Number of rows in the matrix
        n (int): Number of columns in the matrix
        batch_size (int): Number of columns to load and write at a time
        chunk_size (int): Output Zarr column chunk size
    """
    z = zarr.open(zarr_file, mode="a")
    z["X"] = zarr.empty((m, 0), chunks=(m, chunk_size), dtype="float32")

    with h5py.File(file, "r") as f:
        for i in range(n // batch_size + 1):
            j = i * batch_size
            if j + batch_size > n:
                batch_size = n - j
            k = j + batch_size

            matrix = f["X"][:, j:k]

            z["X"].append(matrix, axis=1)
    return


if __name__ == "__main__":
    fire.Fire(h5ad_to_zarr)
