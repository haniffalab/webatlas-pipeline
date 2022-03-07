#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict

import anndata as ad
from shapely import wkt
from shapely.geos import WKTReadingError

from cluster import cluster as get_clusters

def cells_dict(adata):
    cells_dict = {}
    for index, cell in enumerate(adata.obs.index):
        spatial_x, spatial_y, *rest = adata.obsm["spatial"][index]
        try:
            cells_dict[str(adata.obs._indices[index])] = {
                "mappings": {
                    'X_spatial': [float(spatial_x), float(spatial_y)],
                },
                "genes": get_genes(adata, cell),
                "xy": [float(spatial_x), float(spatial_y)],
                "factors": get_factors(adata, index),
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
    genes = [g for g in adata.var_names]
    genes_list = genes[:1000]

    output = {}
    for gene in genes_list:
        try:
            expression = adata[cell,gene].X
            output[gene] = float(expression[0,0])
        except KeyError as e:
            print(e) # @todo
            output[gene] = float(0)
        except IndexError as e:
            print(e) # @todo
            output[gene] = float(0)	

    return output

def get_factors(adata, index):
    factors_list = [
        'sample_name', 
        'total_counts'
    ]
    output = {}

    for factor in factors_list:
        try:
            output[factor] = str(adata.obs[factor][index])
        except KeyError as e:
            print(e) # @todo
            output[factor] = 0
        except IndexError as e:
            print(e) # @todo
            output[factor] = 0	

    return output

def cell_sets_json(data):
    '''
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
    '''
    clustering_dict = {
        'sample_name': defaultdict(list),
        'total_counts': defaultdict(list)
    }
    nice_names = {
        'sample_name': 'Sample Name',
        'total_counts': 'Total Counts'
    }
    for cell_id in data.keys():
        factors = data[cell_id]['factors']
        factors_dict = {
            'sample_name': 'Cluster {}'.format(factors['sample_name'][0]),
            'total_counts': 'Cluster {}'.format(factors['total_counts'][0])
        }
        # For each cluster assignment, append this cell ID to the
        # appropriate clustering_dict list.
        for factor_type, factor_cluster in factors_dict.items():
            clustering_dict[factor_type][factor_cluster].append(cell_id)

    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
    cell_sets = {
        'version': '0.1.2',
        'datatype': 'cell',
        'tree': []
    }

    for factor_type, factor_clusters in clustering_dict.items():
        factor_type_children = []
        for cluster_name in sorted(factor_clusters.keys()):
            factor_type_children.append({
                'name': cluster_name,
                'set': factor_clusters[cluster_name]
            })
        cell_sets['tree'].append({
            'name': nice_names[factor_type],
            'children': factor_type_children
        })

    return cell_sets


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create JSON with cell metadata from H5AD file.')
    parser.add_argument(
        '--h5ad_file', required=True,
        help='H5AD input file.')
    parser.add_argument(
        '--cells_file', type=argparse.FileType('x'),
        help='Write the cell data to this file.')
    parser.add_argument(
        '--cell_sets_file', type=argparse.FileType('x'),
        help='Write the cell-sets data to this file.')
    parser.add_argument(
        '--matrix_file', type=argparse.FileType('x'),
        help='Write the cell-sets data to this file.')
    parser.add_argument(
        '--genes_file', type=argparse.FileType('x'),
        help='Write a list of genes to this file.')
    args = parser.parse_args()

    adata = ad.read(args.h5ad_file)

    metadata = cells_dict(adata)

    if args.cells_file:
        json.dump(metadata, args.cells_file, indent=1)
    if args.cell_sets_file:
        setmetadata = cell_sets_json(metadata)
        json.dump(setmetadata, args.cell_sets_file, indent=1)
    if args.matrix_file:
        clustermetadata = get_clusters(metadata)
        json.dump(clustermetadata, args.matrix_file, indent=1)
