#!/usr/bin/env python3

import json
import argparse
import pickle
from collections import defaultdict

import numpy as np
import pandas

from loom_reader import LoomReader
from cluster import cluster as get_clusters
from delaunay import DictDelaunay2d
from sklearn import decomposition


def octagon(poly):
    '''
    Returns a bounding octagon.

    >>> square = np.array([[0,0], [0,1], [1,1], [1,0]])
    >>> octagon(square)
    [[0, 0], [0, 1], [0, 1], [1, 1], [1, 1], [1, 0], [1, 0], [0, 0]]

    >>> triangle = np.array([[1,0], [0,2], [2,3]])
    >>> octagon(triangle)
    [[0, 1], [0, 2], [1, 3], [2, 3], [2, 3], [2, 1], [1, 0], [1, 0]]

    >>> type(octagon(triangle)[0][0])
    <class 'int'>
    '''
    # SciPy has ConvexHull, but it segfaulted on me: perhaps
    #   https://github.com/scipy/scipy/issues/9751
    # ... and even if I fixed it locally,
    # not sure we want that fragility.
    #
    # Also: The goal is really just to get a simplified shape...
    # A convex hull is too precise in some ways,
    # while in others it falls short, ie, concavities.
    #
    # I kind-of like the obvious artificiality of an octagon.

    # Was unsigned, and substraction causes underflow.
    poly_as_int = poly.astype('int')
    min_x = int(np.min(poly_as_int[:, [0]]))
    max_x = int(np.max(poly_as_int[:, [0]]))
    min_y = int(np.min(poly_as_int[:, [1]]))
    max_y = int(np.max(poly_as_int[:, [1]]))

    summed = np.sum(poly_as_int, axis=1)
    diffed = np.diff(poly_as_int, axis=1)

    min_sum = int(np.min(summed))
    max_sum = int(np.max(summed))
    min_diff = int(np.min(diffed))
    max_diff = int(np.max(diffed))

    return [
        [min_x, min_sum - min_x],
        [min_x, max_diff + min_x],  # ^ Left
        [max_y - max_diff, max_y],
        [max_sum - max_y, max_y],  # ^ Botton
        [max_x, max_sum - max_x],
        [max_x, min_diff + max_x],  # ^ Right
        [min_y - min_diff, min_y],
        [min_sum - min_y, min_y]  # ^ Top
    ]


def mean_coord(coords):
    '''
    The xy values in the Linnarsson data are not good:
    They take a different corner as the origin.
    So... we find the center of our polygon instead.

    >>> mean_coord([[1,2], [3,4], [5,6]])
    [3, 4]

    '''
    return [int(x) for x in np.mean(coords, axis=0).tolist()]


# Taken from http://linnarssonlab.org/osmFISH/clusters/
LOOKUP = {
  "Astrocyte Gfap": "Astrocyte",
  "Astrocyte Mfge8": "Astrocyte",
  "C. Plexus": "Ventricle",
  "Endothelial 1": "Vasculature",
  "Endothelial": "Vasculature",
  "Ependymal": "Ventricle",
  "Hippocampus": "Excitatory neurons",
  "Inhibitory CP": "Inhibitory neurons",
  "Inhibitory Cnr1": "Inhibitory neurons",
  "Inhibitory Crhbp": "Inhibitory neurons",
  "Inhibitory IC": "Inhibitory neurons",
  "Inhibitory Kcnip2": "Inhibitory neurons",
  "Inhibitory Pthlh": "Inhibitory neurons",
  "Inhibitory Vip": "Inhibitory neurons",
  "Microglia": "Brain immune",
  "Oligodendrocyte COP": "Oligodendrocytes",
  "Oligodendrocyte MF": "Oligodendrocytes",
  "Oligodendrocyte Mature": "Oligodendrocytes",
  "Oligodendrocyte NF": "Oligodendrocytes",
  "Oligodendrocyte Precursor cells": "Oligodendrocytes",
  "Pericytes": "Vasculature",
  "Perivascular Macrophages": "Brain immune",
  "Pyramidal Cpne5": "Excitatory neurons",
  "Pyramidal Kcnip2": "Excitatory neurons",
  "Pyramidal L2-3 L5": "Excitatory neurons",
  "Pyramidal L2-3": "Excitatory neurons",
  "Pyramidal L3-4": "Excitatory neurons",
  "Pyramidal L5": "Excitatory neurons",
  "Pyramidal L6": "Excitatory neurons",
  "Vascular Smooth Muscle": "Vasculature",
  "pyramidal L4": "Excitatory neurons"
}


def get_neighborhoods(metadata):
    '''
    >>> cells = {
    ...   'O': { 'xy': [0,0], 'extra': 'field'},
    ...   'N': { 'xy': [0,1], 'extra': 'field'},
    ...   'E': { 'xy': [1,0], 'extra': 'field'},
    ...   'S': { 'xy': [0,-1], 'extra': 'field'},
    ...   'W': { 'xy': [-1,0], 'extra': 'field'}
    ... }
    >>> neighborhoods = get_neighborhoods(cells)
    >>> neighborhoods.keys()
    dict_keys(['O::E::N', 'O::N::W', 'O::S::E', 'O::W::S'])
    >>> neighborhoods['O::E::N']
    {'poly': [[0, 0], [1, 0], [0, 1]]}

    '''
    coords = {}
    for (k, v) in metadata.items():
        coords[k] = v['xy']
    triangles = DictDelaunay2d(coords).getTriangles()
    neighborhoods = {}
    for triangle in triangles:
        key = '::'.join(triangle)
        value = {
            'poly': [coords[point] for point in triangle]
        }
        neighborhoods[key] = value
    return neighborhoods


def get_genes(metadata):
    '''
    >>> metadata = {
    ...   'cell-1': {'genes': {'a': 1, 'b': 20}},
    ...   'cell-2': {'genes': {'a': 2, 'b': 10}}
    ... }
    >>> genes = get_genes(metadata)
    >>> genes['a']
    {'max': 2, 'cells': {'cell-1': 1, 'cell-2': 2}}
    >>> genes['b']
    {'max': 20, 'cells': {'cell-1': 20, 'cell-2': 10}}

    '''
    genes = defaultdict(lambda: {'max': 0, 'cells': {}})
    for cell_id, cell_data in metadata.items():
        for gene_id, expression_level in cell_data['genes'].items():
            gene_data = genes[gene_id]
            gene_data['cells'][cell_id] = expression_level
            if gene_data['max'] < expression_level:
                gene_data['max'] = expression_level
    return genes


def get_factors(metadata):
    '''
    >>> metadata = {
    ...   "Santa's Little Helper": {'factors':{'eng': 'dog', 'sci': 'canine'}},
    ...   "Snowball II": {'factors':{'eng': 'cat', 'sci': 'feline'}}
    ... }
    >>> factors = get_factors(metadata)
    >>> list(factors['eng'].keys())
    ['map', 'cells']
    >>> factors['eng']['map']
    ['dog', 'cat']
    >>> factors['eng']['cells']
    {"Santa's Little Helper": 0, 'Snowball II': 1}

    '''
    factors = defaultdict(lambda: {'map': [], 'cells': {}})
    for cell_id, cell_data in metadata.items():
        for factor_id, factor_value in cell_data['factors'].items():
            factor_data = factors[factor_id]
            if factor_value not in factor_data['map']:
                factor_data['map'].append(factor_value)
            factor_index = factor_data['map'].index(factor_value)
            factor_data['cells'][cell_id] = factor_index
    return factors


def get_cell_sets(clusters, lookup):
    '''
    >>> from collections import namedtuple
    >>> Cluster = namedtuple('Cluster', ['name', 'cell_ids'])
    >>> clusters = {
    ...   1: Cluster('pyramidal L4', ['104', '110', '111']),
    ...   3: Cluster('vascular smooth muscle', ['1', '2', '3'])
    ... }
    >>> lookup = {
    ...   'vascular smooth muscle': 'vasculature',
    ...   'pyramidal L4': 'excitatory neurons'
    ... }
    >>> cell_sets = get_cell_sets(clusters, lookup)
    >>> list(cell_sets.keys())
    ['version', 'datatype', 'tree']
    >>> cell_sets['datatype']
    'cell'
    >>> list(cell_sets['tree'][0].keys())
    ['name', 'children']
    >>> cell_sets['tree'][0]['name']
    'Cell Type Annotations'
    >>> sorted([ n['name'] for n in cell_sets['tree'][0]['children'] ])
    ['excitatory neurons', 'vasculature']
    '''

    # The parameter `lookup` is a dict mapping
    # subclusters to clusters: `{ Subcluster Name: Cluster Name }`
    # This `lookup` mapping can be used to fill in an intermediate
    # dict `hierarchy`, closer to the data structure we want to output.
    # ```
    # {
    #   Cluster A: {
    #     Subcluster A: [1, 2],
    #     Subcluster B: [3, 4]
    #   },
    #   Cluster B: {...}
    # }
    # ```
    hierarchy = {cluster_name: {} for cluster_name in lookup.values()}
    for c in clusters.values():
        subcluster = {c.name: c.cell_ids}
        cluster_name = lookup.get(c.name)
        cluster_dict = hierarchy.get(cluster_name)
        cluster_dict.update(subcluster)

    # Use the `hierarchy` dict to fill in an object
    # conforming to the `cell-sets.json` schema.
    cluster_nodes = []
    for cluster_name in sorted(hierarchy.keys()):
        cluster_dict = hierarchy[cluster_name]
        subcluster_nodes = []
        for subcluster_name in sorted(cluster_dict.keys()):
            subcluster = cluster_dict[subcluster_name]
            subcluster_nodes.append({
                'name': subcluster_name,
                'set': subcluster
            })
        cluster_nodes.append({
            'name': cluster_name,
            'children': subcluster_nodes,
        })

    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
    cell_sets = {
        'version': '0.1.2',
        'datatype': 'cell',
        'tree': [{
            'name': 'Cell Type Annotations',
            'children': cluster_nodes
        }]
    }

    return cell_sets


def genes_to_samples_by_features(metadata):
    '''
    >>> metadata = {
    ...   '0': {
    ...     'genes': {'A': 0, 'B': 0, 'A2': 0, 'B2': 0}
    ...   },
    ...   '1': {
    ...     'genes': {'A': 0, 'B': 1, 'A2': 0, 'B2': 1}
    ...   },
    ...   '2': {
    ...     'genes': {'A': 0, 'B': 4, 'A2': 0, 'B2': 4}
    ...   }
    ... }
    >>> s_by_f = genes_to_samples_by_features(metadata)
    >>> s_by_f.shape
    (3, 4)
    '''
    records = dict([(k, v['genes']) for k, v in metadata.items()])
    return pandas.DataFrame.from_dict(records, orient='index')


def add_pca(metadata):
    '''
    >>> metadata = {
    ...   '0': {
    ...     'mappings': {},
    ...     'genes': {'A': 0, 'B': 0, 'A2': 1, 'B2': 0}
    ...   },
    ...   '1': {
    ...     'mappings': {},
    ...     'genes': {'A': 1, 'B': 1, 'A2': 0, 'B2': 1}
    ...   },
    ...   '2': {
    ...     'mappings': {},
    ...     'genes': {'A': 0, 'B': 4, 'A2': 0, 'B2': 4}
    ...   }
    ... }
    >>> add_pca(metadata)
    >>> metadata['0']['mappings']['PCA']
    [-2.41, -0.57]
    >>> metadata['1']['mappings']['PCA']
    [-0.92, 0.77]
    >>> metadata['2']['mappings']['PCA']
    [3.33, -0.2]
    '''
    pca = decomposition.PCA(n_components=2)
    principle_components = pca.fit_transform(
        genes_to_samples_by_features(metadata)
    ).tolist()
    for (k, pc) in zip(metadata.keys(), principle_components):
        metadata[k]['mappings']['PCA'] = [
            round(component, 2) for component in pc
        ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create JSON with cell metadata and, '
                    'optionally, segmentation.')
    parser.add_argument(
        '--loom', required=True,
        help='Loom file with cell metadata')
    parser.add_argument(
        '--pkl', type=argparse.FileType('rb'),
        help='Pickle file with cell segmentation data')
    parser.add_argument(
        '--clusters_file', type=argparse.FileType('x'),
        help='Write the hierarchically clustered data to this file.')
    parser.add_argument(
        '--cells_file', type=argparse.FileType('x'),
        help='Write the cleaned cell data to this file.')
    parser.add_argument(
        '--cell-sets_file', type=argparse.FileType('x'),
        help='Write the cleaned cell sets data to this file.')
    parser.add_argument(
        '--genes_file', type=argparse.FileType('x'),
        help='Write a list of genes to this file.'),
    parser.add_argument(
        '--neighborhoods_file', type=argparse.FileType('x'),
        help='Write the cell neighborhoods to this file.')
    parser.add_argument(
        '--factors_file', type=argparse.FileType('x'),
        help='Write the cell factors to this file.')
    parser.add_argument(
        '--integers', action='store_true',
        help='Convert all numbers to integers.')
    args = parser.parse_args()

    lr = LoomReader(args.loom)
    metadata = lr.data()
    add_pca(metadata)

    for cell in metadata.values():
        # "Clusters" in the raw data are called "subclusters"
        # in http://linnarssonlab.org/osmFISH/clusters/
        subcluster = cell.pop('cluster')
        cell['factors'] = {
            'subcluster': subcluster,
            'cluster': LOOKUP[subcluster]
        }

    if args.pkl:
        segmentation = pickle.load(args.pkl)
        for cell_id, poly in segmentation.items():
            if cell_id in metadata:
                simple_poly = octagon(poly)
                xy = mean_coord(simple_poly)
                metadata[cell_id]['poly'] = simple_poly
                metadata[cell_id]['xy'] = xy

    if args.integers:
        for cell in metadata.values():
            cell['xy'] = [
                # Raw data has way too many decimal points!
                int(z) for z in cell['xy']
            ]

    if args.cells_file:
        json.dump(metadata, args.cells_file, indent=1)

    if args.cell_sets_file:
        clusters = lr.clusters()
        cell_sets = get_cell_sets(clusters, LOOKUP)
        json.dump(cell_sets, args.cell_sets_file, indent=1)

    if args.clusters_file:
        clusters = get_clusters(metadata)
        clusters_json = json.dumps(clusters)
        # Line-break after every element is too much, but this works:
        spaced_clusters_json = clusters_json.replace(
            '],',
            '],\n'
        )
        print(spaced_clusters_json, file=args.clusters_file)

    if args.genes_file:
        genes = get_genes(metadata)
        genes_json = json.dumps(genes)
        spaced_genes_json = genes_json.replace(
            '},',
            '},\n'
        )
        print(spaced_genes_json, file=args.genes_file)

    if args.factors_file:
        factors = get_factors(metadata)
        factors_json = json.dumps(factors)
        spaced_factors_json = factors_json.replace(
            '},',
            '},\n'
        )
        print(spaced_factors_json, file=args.factors_file)

    if args.neighborhoods_file:
        neighborhoods = get_neighborhoods(metadata)
        json.dump(neighborhoods, args.neighborhoods_file)
