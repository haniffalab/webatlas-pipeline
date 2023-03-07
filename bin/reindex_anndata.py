import os
import zarr
import anndata as ad
from scipy.sparse import spmatrix


def reindex_anndata(
    anndata_file: str,
    offset: int,
    chunk_size: int = 10,
    override_chunk_size: bool = False,
):

    basename, ext = os.path.splitext(anndata_file)
    out_filename = f"reindexed-{basename}"
    is_zarr = ext == ".zarr"

    if is_zarr:
        z = zarr.open(anndata_file)
        adata = ad.read_zarr(z.store)
    else:
        adata = ad.read(anndata_file)

    adata.obs.index = (adata.obs.index.astype(int) + offset).astype(str)

    if is_zarr and not override_chunk_size:
        chunk_shape = z.X.chunks
    else:
        chunk_shape = (adata.shape[0], chunk_size)

    if isinstance(adata.X, spmatrix):
        adata.X = adata.X.toarray()

    adata.write_zarr(f"{out_filename}.zarr", chunk_shape)
