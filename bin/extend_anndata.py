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
    stem = stem or "{}-extended".format(basename)
    is_zarr = ext == ".zarr"
    if is_zarr:
        z = zarr.open(anndata_file)
        adata = ad.read_zarr(z.store)
    else:
        adata = ad.read(anndata_file)

    ext_matrix = pd.get_dummies(adata.obs[obs], dtype="float32")

    obs_feature_name = obs_feature_name or obs
    prev_features_bool = "is_{}".format(feature_name)
    new_features_bool = "is_{}".format(obs_feature_name)

    if isinstance(adata.X, spmatrix):
        adata_combined = ad.AnnData(
            hstack(
                (
                    adata.X,
                    csr_matrix(ext_matrix.values)
                    if isinstance(adata.X, csr_matrix)
                    else csc_matrix(ext_matrix.values),
                )
            ),
            obs=adata.obs,
            var=pd.concat(
                [
                    adata.var.assign(**{prev_features_bool: True}),
                    ext_matrix.columns.to_frame(obs_feature_name)
                    .drop(columns=0)
                    .assign(**{new_features_bool: True}),
                ]
            ),
            obsm=adata.obsm,
            uns=adata.uns,
        )
    else:
        adata_combined = ad.AnnData(
            np.hstack((adata.X, ext_matrix.values)),
            obs=adata.obs,
            var=pd.concat(
                [
                    adata.var.assign(**{prev_features_bool: True}),
                    ext_matrix.columns.to_frame(obs_feature_name)
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
    
    if is_zarr and not override_chunk_size:
        chunk_shape = z.X.chunks
    else:
        chunk_shape = (adata.shape[0], chunk_size)
    
    zarr_file = f"{stem}-{SUFFIX}"

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
    stem = stem or "{}-extended".format(basename)
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

    cell_types = ad.AnnData(
        c2l_adata.obsm[q].to_numpy(),
        dtype="float32",
        obs=c2l_adata.obs,
        var=c2l_adata.obsm[q].columns.str.replace(q.split("_")[0] + "cell_abundance_w_sf_", "").to_frame(),
    )
    cell_types.var.drop(columns=0, inplace=True)
    cell_types.var.index.rename("CellType", inplace=True)

    if sample:
        cell_types = cell_types[cell_types.obs[sample[0]] == sample[1]]

    
    # Create AnnData concatenating gene expression matrix and cell type abundance matrix
    if isinstance(adata.X, spmatrix):
        adata_combined = ad.AnnData(
            hstack((
                    adata.X,
                    csr_matrix(cell_types.X)
                    if isinstance(adata.X, csr_matrix)
                    else csc_matrix(cell_types.X),
            )),
            obs=adata.obs,
            var=pd.concat([
                    adata.var.assign(is_gene=True),
                    cell_types.var.assign(is_celltype=True),
                ],
                axis=1,
            ),
        )
    else:
        adata_combined = ad.AnnData(
            np.hstack((adata.X, cell_types.X)),
            obs=adata.obs,
            var=pd.concat([
                    adata.var.assign(is_gene=True),
                    cell_types.var.assign(is_celltype=True),
                ],
                axis=1,
            ),
        )

    adata_combined.var["is_gene"] = adata_combined.var["is_gene"].fillna(False)
    adata_combined.var["is_celltype"] = adata_combined.var["is_celltype"].fillna(False)


    # Ensure bool columns remain bool, filling nans with False or fill as appropriate
    for col in [col for col in adata.var_keys() if adata.var[col].dtype == bool]:
        adata_combined.var[col] = adata_combined.var[col].fillna(False)
    for col in [col for col in cell_types.var_keys() if cell_types.var[col].dtype == bool]:
        adata_combined.var[col] = adata_combined.var[col].fillna(False)
    
    if is_zarr and not override_chunk_size:
        chunk_shape = z.X.chunks
    else:
        chunk_shape = (adata.shape[0], chunk_size)
    
    zarr_file = f"{stem}-{SUFFIX}"

    if isinstance(adata.X, spmatrix):
        adata.X = adata.X.toarray()

    adata.write_zarr(zarr_file, chunk_shape)
