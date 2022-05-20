from vitessce import (
    DataType as dt,
    FileType as ft,
    Component as cm,
    CoordinationType as ct
)

SINGLE_ZARR = 'anndata.zarr'

# Data types with ordered file types
DATA_TYPES = {
    dt.CELLS: [
        ('cells.json', ft.CELLS_JSON),
        ('anndata-cells.zarr', ft.ANNDATA_CELLS_ZARR),
        (SINGLE_ZARR, ft.ANNDATA_CELLS_ZARR),
    ],
    dt.MOLECULES: [
        ('molecules.json', ft.MOLECULES_JSON),
    ],
    dt.CELL_SETS: [
        ('cell-sets.json', ft.CELL_SETS_JSON),
        ('anndata-cell-sets.zarr', ft.ANNDATA_CELL_SETS_ZARR),
        (SINGLE_ZARR, ft.ANNDATA_CELL_SETS_ZARR),
    ],
    dt.RASTER: [
        ('raster.ome-zarr', 'raster.ome-zarr'),
        ('raster.json', ft.RASTER_JSON),
    ],
    dt.EXPRESSION_MATRIX: [
        ('expression-matrix.zarr', ft.EXPRESSION_MATRIX_ZARR),
        ('anndata-expression-matrix.zarr', ft.ANNDATA_EXPRESSION_MATRIX_ZARR),
        ('clusters.json', ft.CLUSTERS_JSON),
        ('genes.json', ft.GENES_JSON),
        (SINGLE_ZARR, ft.ANNDATA_EXPRESSION_MATRIX_ZARR),
    ],
    dt.NEIGHBORHOODS: [
        ('neighborhoods.json', ft.NEIGHBORHOODS_JSON),
    ],
    dt.GENOMIC_PROFILES: [
        ('genomic-profiles.zarr', ft.GENOMIC_PROFILES_ZARR),
    ],
}

DEFAULT_OPTIONS = {
    ft.ANNDATA_CELLS_ZARR: {
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
    ft.ANNDATA_CELL_SETS_ZARR: {
        'sets': [ 'obs/sample' ]
    },
    ft.ANNDATA_EXPRESSION_MATRIX_ZARR: {
        'matrix': 'X'
    }
}

DEFAULT_LAYOUTS = {
    'minimal': [
        [cm.SPATIAL],
        [cm.LAYER_CONTROLLER]
    ],
    'simple': [
        [cm.SPATIAL],
        [cm.LAYER_CONTROLLER],
        [cm.GENES],
        [cm.CELL_SETS]
    ],
    'advanced': [
        [cm.LAYER_CONTROLLER],
        [cm.SPATIAL],
        [cm.SCATTERPLOT, cm.SCATTERPLOT],
        [cm.GENES, cm.CELL_SETS],
        [cm.GENOMIC_PROFILES]
    ]
}

# Coordination Types required by Components/Views
COMPONENTS_COORDINATION_TYPES = {
    cm.SCATTERPLOT : [
        ct.EMBEDDING_TYPE
    ]
}

# Data Types required by Components/Views
COMPONENTS_DATA_TYPES = {
    cm.SCATTERPLOT: set([
        dt.CELLS
    ]),
    cm.HEATMAP: set([
        dt.EXPRESSION_MATRIX
    ]),
    cm.SPATIAL: set([
        dt.RASTER,
        dt.CELLS,
        dt.MOLECULES
    ]),
    cm.LAYER_CONTROLLER: set([
        dt.RASTER,
        dt.CELLS,
        dt.MOLECULES
    ]),
    cm.GENOMIC_PROFILES: set([
        dt.GENOMIC_PROFILES
    ]),
    cm.GENES: set([
        dt.EXPRESSION_MATRIX
    ]),
    cm.CELL_SETS: set([
        dt.CELL_SETS
    ]),
    cm.CELL_SET_SIZES: set([
        dt.CELL_SETS
    ]),
    cm.CELL_SET_EXPRESSION: set([
        dt.EXPRESSION_MATRIX
    ])
}