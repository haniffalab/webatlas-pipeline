#!/usr/bin/env python3

# import argparse
import fire
import json
from collections import defaultdict

import anndata as ad

# from shapely import wkt
# from shapely.geos import WKTReadingError

from cluster import cluster as get_clusters


def cells_dict(adata):
    cells_dict = {}
    for index, cell in enumerate(adata.obs.index):
        # polygon = adata.obsm["polys"][index]
        # simpler = adata.obsm["simpler"][index]
        spatial = adata.obsm["spatial"][index]
        # pca_x, pca_y, *rest = adata.obsm["X_pca"][index]
        # umap_x, umap_y = adata.obsm["X_umap"][index]
        try:
            # poly = list(wkt.loads(simpler).coords)
            cells_dict[cell] = {
                "mappings": {
                    # "X_pca": [float(pca_x), float(pca_y)],
                    # "X_umap": [float(umap_x), float(umap_y)],
                },
                "genes": get_genes(adata, cell),
                "xy": list(spatial),
                "factors": get_factors(adata, index),
                # "poly": poly,
            }
        except ValueError as e:
            print(e)
            pass
        except WKTReadingError as e:
            print(e)
            pass
        except NotImplementedError as e:
            print(e)
            pass

    return cells_dict


def get_genes(adata, cell):
    genes_list = adata.var.index[:10]
    output = {}

    for gene in genes_list:
        try:
            expression = adata[cell, gene].X
            output[gene] = float(expression[0, 0])
        except KeyError as e:
            print(e)  # @todo
            output[gene] = float(0)
        except IndexError as e:
            print(e)  # @todo
            output[gene] = float(0)

    return output


def get_factors(adata, index):
    factors_list = ["sample", "total_counts", "n_genes_by_counts", "doublet_scores"]
    output = {}

    for factor in factors_list:
        try:
            output[factor] = str(adata.obs[factor][index])
        except KeyError as e:
            print(e)  # @todo
            output[factor] = 0
        except IndexError as e:
            print(e)  # @todo
            output[factor] = 0

    return output


def cell_sets_json(data):
    """
    >>> data = {
    ...     "cell_1": {
    ...         "locations": [1632.02, -1305.7],
    ...         "factors": {
    ...             "pleiden_clus": [3],
    ...             "kmeans": [8]
    ...         },
    ...         "mappings": {
    ...             "tsne": [-11.0776, 6.0311],
    ...             "umap": [-6.858, 15.7691]
    ...         }
    ...     }
    ... }
    >>> cell_sets = cell_sets_json(data)
    >>> list(cell_sets.keys())
    ['version', 'datatype', 'tree']
    >>> cell_sets['datatype']
    'cell'
    >>> len(cell_sets['tree'])
    2
    >>> sorted([ n['name'] for n in cell_sets['tree'] ])
    ['Leiden Clustering', 'k-means Clustering']
    """

    clustering_dict = {"sample": defaultdict(list), "total_counts": defaultdict(list)}
    nice_names = {"sample": "Sample", "total_counts": "Total Counts"}
    for cell_id in data.keys():
        factors = data[cell_id]["factors"]
        factors_dict = {
            "sample": "Cluster {}".format(factors["sample"][0]),
            "total_counts": "Cluster {}".format(factors["total_counts"][0]),
        }
        # For each cluster assignment, append this cell ID to the
        # appropriate clustering_dict list.
        for factor_type, factor_cluster in factors_dict.items():
            clustering_dict[factor_type][factor_cluster].append(cell_id)

    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
    cell_sets = {"version": "0.1.2", "datatype": "cell", "tree": []}

    for factor_type, factor_clusters in clustering_dict.items():
        factor_type_children = []
        for cluster_name in sorted(factor_clusters.keys()):
            factor_type_children.append(
                {"name": cluster_name, "set": factor_clusters[cluster_name]}
            )
        cell_sets["tree"].append(
            {"name": nice_names[factor_type], "children": factor_type_children}
        )

    return cell_sets


def main(h5ad_file, cells_file=None, cell_sets_file=None, matrix_file=None):
    metadata = cells_dict(ad.read(h5ad_file))
    # print(metadata)

    if cells_file:
        json.dump(metadata, open(cells_file, "w"), indent=1)
    if cell_sets_file:
        json.dump(cell_sets_json(metadata), open(cell_sets_file, "w"), indent=1)
    if matrix_file:
        json.dump(get_clusters(metadata), open(matrix_file, "w"), indent=1)


if __name__ == "__main__":
    fire.Fire(main)
