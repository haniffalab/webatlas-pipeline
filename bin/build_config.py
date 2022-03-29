#!/usr/bin/env python3

import os
import fire
import json

DATA_TYPES = {
    'cells': [
        'cells.json',
        'anndata-cells.zarr',
    ],
    'molecules': [
        'molecules.json',
    ],
    'cell-sets': [
        'cell-sets.json',
        'anndata-cell-sets.zarr',
    ],
    'raster': [
        'raster.ome-zarr',
        'raster.json',
    ],
    'expression-matrix': [
        'expression-matrix.zarr',
        'anndata-expression-matrix.zarr',
        'clusters.json',
        'genes.json',
    ],
    'neighborhoods': [
        'neighborhoods.json'
    ],
    'genomic-profiles': [
        'genomic-profiles.zarr'
    ],
}

def write_raster_json(
    dataset='',
    files_dir='',
    zarr_dirs=None,
    url='',
    outdir='',
    config_filename='config.json',
    ):
    config_json = {
        "name": dataset,
        "version": "1.0.0",
        "description": "",
        "public": True,
        "datasets": [],
        "initStrategy": "auto",
        "coordinationSpace": {},
        "layout": []
    }

    files = []
    for data_type in DATA_TYPES:
        for file_type in DATA_TYPES[data_type]:
            f = os.path.join(files_dir, file_type)
            if os.path.isfile(f):
                files.append({
                    "url": file_type,
                    "type": data_type,
                    "fileType": file_type
                })

    # TODO: add zarr data
    if zarr_dirs and type(zarr_dirs) != bool:
        if type(zarr_dirs) not in [list, tuple]:
            zarr_dirs = [zarr_dirs]
    
    config_json["datasets"].append({
        "uid": dataset,
        "name": dataset,
        "description": "",
        "files": files,
    })

    if outdir and not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir or '', config_filename), "w") as out_file:
        json.dump(config_json, out_file, indent=2)


if __name__ == "__main__":
    fire.Fire(write_raster_json)