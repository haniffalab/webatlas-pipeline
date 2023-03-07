import os
import zarr
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import spmatrix, hstack, csr_matrix, csc_matrix


SUFFIX = "anndata.zarr"

def from_obs(
    anndata_file: str,
    obs: str = "celltype",
    stem: str = None,
    feature_name: str = "gene",
    obs_feature_name: str = None,
    chunk_size: int = 10,
    override_chunk_size: bool = False
):

    basename, ext = os.path.splitext(anndata_file)
    out_filename = f"{stem}-{SUFFIX}" if stem else f"{basename}-extended"
    zarr_file = out_filename + ".zarr"
    is_zarr = ext == ".zarr"

    if is_zarr:
        z = zarr.open(anndata_file)
        adata = ad.read_zarr(z.store)
    else:
        adata = ad.read(anndata_file)

    ext_matrix = pd.get_dummies(adata.obs[obs], dtype="float32")

    adata = extend_matrix(adata, ext_matrix)
    
    if is_zarr and not override_chunk_size:
        chunk_shape = z.X.chunks
    else:
        chunk_shape = (adata.shape[0], chunk_size)

    if isinstance(adata.X, spmatrix):
        adata.X = adata.X.toarray()

    adata.write_zarr(zarr_file, chunk_shape)


def from_cell2location(
    anndata_file: str,
    c2l_file: str,
    stem: str = None,
    q: str = "q05_cell_abundance_w_sf",
    sample: str = None,
    chunk_size: int = 10,
    override_chunk_size: bool = False
):

    basename, ext = os.path.splitext(anndata_file)
    out_filename = f"{stem}-{SUFFIX}" if stem else f"{basename}-extended"
    zarr_file = out_filename + ".zarr"
    is_zarr = ext == ".zarr"
    
    if is_zarr:
        z = zarr.open(anndata_file)
        adata = ad.read_zarr(z.store)
    else:
        adata = ad.read(anndata_file)
    

    # Read only obs, var and obsm from the cell2location output h5ad
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
        columns=c2l_adata.obsm[q].columns.str.replace(q.split("_")[0] + "cell_abundance_w_sf_", ""),
        dtype="float32"
    )

    adata = extend_matrix(adata, c2l_df)
    
    if is_zarr and not override_chunk_size:
        chunk_shape = z.X.chunks
    else:
        chunk_shape = (adata.shape[0], chunk_size)

    if isinstance(adata.X, spmatrix):
        adata.X = adata.X.toarray()

    adata.write_zarr(zarr_file, chunk_shape)


def extend_matrix(
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
        adata_combined = ad.AnnData(
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
        adata_combined = ad.AnnData(
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

    adata_combined.var[prev_features_bool] = adata_combined.var[
        prev_features_bool
    ].fillna(False)
    adata_combined.var[new_features_bool] = adata_combined.var[
        new_features_bool
    ].fillna(False)

    for col in [col for col in adata.var_keys() if adata.var[col].dtype == bool]:
        adata_combined.var[col] = adata_combined.var[col].fillna(False)
    
    return adata_combined