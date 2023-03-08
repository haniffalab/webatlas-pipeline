#!/usr/bin/env python3

from typing import Union
import os
import fire
import zarr
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import spmatrix, hstack, csr_matrix, csc_matrix
from process_h5ad import h5ad_to_zarr


def reindex_anndata(path: str, offset: int, **kwargs):

    adata = read_anndata(path)

    adata.obs.index = (adata.obs.index.astype(int) + offset).astype(str)

    out_filename = "reindexed-{}".format(os.path.splitext(os.path.basename(path))[0])

    write_anndata(adata, out_filename, **kwargs)

    return


def concat_features(path: str, features: str, **kwargs):

    adata = read_anndata(path)

    if features.endswith(".h5ad") and os.path.isfile(features):
        adata = concat_matrix_from_cell2location(adata, features)
    else:
        adata = concat_matrix_from_obs(adata, features)

    out_filename = "concat-{}".format(os.path.splitext(os.path.basename(path))[0])

    write_anndata(adata, out_filename, **kwargs)


def intersect_features(paths: list[str], **kwargs):

    var_intersect = get_feature_intersection(paths)

    for path in paths:
        adata = read_anndata(path)

        adata = adata[:, var_intersect]

        out_filename = "intersect-{}".format(
            os.path.splitext(os.path.basename(path))[0]
        )

        write_anndata(adata, out_filename, **kwargs)

    return


def concat_matrix_from_obs(
    data: Union[ad.AnnData, str],
    obs: str = "celltype",
    feature_name: str = "gene",
    obs_feature_name: str = None,
):
    if isinstance(data, ad.AnnData):
        adata = data
    else:
        adata = read_anndata(data)

    ext_matrix = pd.get_dummies(adata.obs[obs], dtype="float32")

    return concat_matrices(adata, ext_matrix, obs, feature_name, obs_feature_name)


def concat_matrix_from_cell2location(
    data: Union[ad.AnnData, str],
    c2l_file: str,
    q: str = "q05_cell_abundance_w_sf",
    sample: str = None,
    feature_name: str = "gene",
    obs_feature_name: str = None,
):

    if isinstance(data, ad.AnnData):
        adata = data
    else:
        adata = read_anndata(data)

    with h5py.File(c2l_file) as f:
        c2l_adata = ad.AnnData(
            obs=ad._io.h5ad.read_elem(f["obs"]) if "obs" in f else None,
            var=ad._io.h5ad.read_elem(f["var"]) if "var" in f else None,
            obsm=ad._io.h5ad.read_elem(f["obsm"]) if "obsm" in f else None,
        )

    if sample:
        c2l_adata = c2l_adata[c2l_adata.obs[sample[0]] == sample[1]]

    c2l_df = pd.DataFrame(
        c2l_adata.obsm[q].to_numpy(),
        index=c2l_adata.obs.index,
        columns=c2l_adata.obsm[q].columns.str.replace(
            q.split("_")[0] + "cell_abundance_w_sf_", ""
        ),
        dtype="float32",
    )

    return concat_matrices(adata, c2l_df, "celltype", feature_name, obs_feature_name)


def concat_matrices(
    adata: ad.AnnData,
    ext_df: pd.DataFrame,
    obs: str = "celltype",
    feature_name: str = "gene",
    obs_feature_name: str = None,
):

    assert adata.shape[0] == ext_df.shape[0]

    obs_feature_name = obs_feature_name or obs
    prev_features_bool = "is_{}".format(feature_name)
    new_features_bool = "is_{}".format(obs_feature_name)

    if isinstance(adata.X, spmatrix):
        adata_concat = ad.AnnData(
            hstack(
                (
                    adata.X,
                    csr_matrix(ext_df.values)
                    if isinstance(adata.X, csr_matrix)
                    else csc_matrix(ext_df.values),
                )
            ),
            obs=adata.obs,
            var=pd.concat(
                [
                    adata.var.assign(**{prev_features_bool: True}),
                    ext_df.columns.to_frame(obs_feature_name)
                    .drop(columns=0)
                    .assign(**{new_features_bool: True}),
                ]
            ),
            obsm=adata.obsm,
            uns=adata.uns,
        )
    else:
        adata_concat = ad.AnnData(
            np.hstack((adata.X, ext_df.values)),
            obs=adata.obs,
            var=pd.concat(
                [
                    adata.var.assign(**{prev_features_bool: True}),
                    ext_df.columns.to_frame(obs_feature_name)
                    .drop(columns=0)
                    .assign(**{new_features_bool: True}),
                ]
            ),
            obsm=adata.obsm,
            uns=adata.uns,
        )

    adata_concat.var[prev_features_bool] = adata_concat.var[prev_features_bool].fillna(
        False
    )
    adata_concat.var[new_features_bool] = adata_concat.var[new_features_bool].fillna(
        False
    )

    for col in [col for col in adata.var_keys() if adata.var[col].dtype == bool]:
        adata_concat.var[col] = adata_concat.var[col].fillna(False)

    return adata_concat


def get_feature_intersection(paths: list[str]):

    var_indices = []
    for path in paths:
        is_zarr = os.path.splitext(path)[-1] == ".zarr"
        if is_zarr:
            z = zarr.open(path)
            var_indices.append(pd.Index(z.var._index[:]).to_series())
        else:
            with h5py.File(path, "r") as f:
                var_indices.append(ad._io.h5ad.read_elem(f["var"]).index.to_series())

    var_intersect = pd.concat(var_indices, axis=1, join="inner").index

    return var_intersect


def read_anndata(path: str):
    is_zarr = os.path.splitext(path)[-1] == ".zarr"

    if is_zarr:
        z = zarr.open(path)
        adata = ad.read_zarr(z.store)
    else:
        adata = ad.read(path)

    return adata


def write_anndata(
    adata: ad.AnnData, out_filename: str, save_h5ad: bool = False, **kwargs
):
    if save_h5ad:
        adata.write_h5ad(f"{out_filename}.h5ad")

    h5ad_to_zarr(adata=adata, out_filename=out_filename, **kwargs)

    return


if __name__ == "__main__":
    fire.Fire()
