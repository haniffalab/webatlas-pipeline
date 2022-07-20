#!/usr/bin/env python3

import fire
import os

from process_h5ad import h5ad_to_zarr
from process_molecules import tsv_to_json
import scanpy as sc
import pandas as pd


def load_space_ranger(dir_path):
    adata = sc.read_visium(dir_path)
    graphclust = pd.read_csv(f"{dir_path}/analysis/clustering/graphclust/clusters.csv", index_col="Barcode")
    adata.obs["graphclust"] = graphclust.Cluster.astype("category")
    umap = pd.read_csv(f"{dir_path}/analysis/umap/2_components/projection.csv", index_col="Barcode")
    tsne = pd.read_csv(f"{dir_path}/analysis/tsne/2_components/projection.csv", index_col="Barcode")
    pca = pd.read_csv(f"{dir_path}/analysis/pca/10_components/projection.csv", index_col="Barcode")
    adata.obsm["X_umap"] = umap.values
    adata.obsm["X_tsne"] = tsne.values
    adata.obsm["X_pca"] = pca.values
    return adata


def main(file_type, file_path, stem, args):
  out_file = None

  if os.path.isdir(file_path):
      adata = load_space_ranger(file_path)
      # print(adata)
      out_file = h5ad_to_zarr(adata=adata, stem=stem, **args)
  elif file_type == 'h5ad':
      out_file = h5ad_to_zarr(adata=sc.read_h5ad(file_path), stem=stem, **args)
  elif file_type == 'molecules':
      out_file = tsv_to_json(file=file_path, stem=stem, **args)

  return out_file

if __name__ == "__main__":
    fire.Fire(main)
