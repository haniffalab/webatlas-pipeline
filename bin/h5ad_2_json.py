#!/usr/bin/env python3

import anndata as ad
from collections import defaultdict
import fire
import json
import pickle

from cluster import cluster as get_clusters


def cells_dict(adata):
    cells_dict = {}
    for index, cell in enumerate(adata.obs.index):
        mappings = {}
        for i in adata.obsm:
            try:
                x, y, *rest = adata.obsm[i][index]
                mappings[i] = [float(x),float(y)]
            except ValueError as e:
                print(e) # @todo
                pass
        try:
            cells_dict[cell] = {
                "mappings": mappings,
                "genes": get_genes(adata, cell),
                "xy": list(map(int, mappings["spatial"]) if "spatial" in mappings else [0,0]),
                "factors": get_factors(adata, index),
            }
        except ValueError as e:
            print(e) # @todo
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
            print(e) # @todo
            output[gene] = float(0)
        except IndexError as e:
            print(e) # @todo
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
    with open('metadata.pickle', 'wb') as handle:
        pickle.dump(metadata, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if cells_file:
        json.dump(metadata, open(cells_file, "w"), indent=1)
    if cell_sets_file:
        json.dump(cell_sets_json(metadata), open(cell_sets_file, "w"), indent=1)
    if matrix_file:
        json.dump(get_clusters(metadata), open(matrix_file, "w"), indent=1)


if __name__ == "__main__":
   fire.Fire(main)
