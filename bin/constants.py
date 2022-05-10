SINGLE_ZARR = 'anndata.zarr'

# Data types with ordered file types
DATA_TYPES = {
    'cells': [
        ('cells.json', 'cells.json'),
        ('anndata-cells.zarr', 'anndata-cells.zarr'),
        (SINGLE_ZARR, 'anndata-cells.zarr'),
    ],
    'molecules': [
        ('molecules.json', 'molecules.json'),
    ],
    'cell-sets': [
        ('cell-sets.json', 'cell-sets.json'),
        ('anndata-cell-sets.zarr', 'anndata-cell-sets.zarr'),
        (SINGLE_ZARR, 'anndata-cell-sets.zarr'),
    ],
    'raster': [
        ('raster.ome-zarr', 'raster.ome-zarr'),
        ('raster.json', 'raster.json'),
    ],
    'expression-matrix': [
        ('expression-matrix.zarr', 'expression-matrix.zarr'),
        ('anndata-expression-matrix.zarr', 'anndata-expression-matrix.zarr'),
        ('clusters.json', 'clusters.json'),
        ('genes.json', 'genes.json'),
        (SINGLE_ZARR, 'anndata-expression-matrix.zarr'),
    ],
    'neighborhoods': [
        ('neighborhoods.json', 'neighborhoods.json'),
    ],
    'genomic-profiles': [
        ('genomic-profiles.zarr', 'genomic-profiles.zarr'),
    ],
}

DEFAULT_OPTIONS = {
    'anndata-cells.zarr': {
        'spatial': {
            'xy': 'obsm/spatial',
        },
        'mappings': {
            'spatial' : None,
        },
        'factors': [
            'obs/sample',
        ]
    },
    'anndata-cell-sets.zarr': {
        'sets': [ 'obs/sample' ]
    },
    'anndata-expression-matrix.zarr': {
        'matrix': 'X'
    }
}
