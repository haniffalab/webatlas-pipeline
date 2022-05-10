#!/usr/bin/env python3

import os
import fire
import json

from constants import DATA_TYPES, DEFAULT_OPTIONS

def build_options(file_type, file_path, file_options=None, check_exist=True):
    options = None
    if file_options is None and file_type in DEFAULT_OPTIONS:
        file_options = DEFAULT_OPTIONS[file_type]
    
    if file_type == 'anndata-cells.zarr':
        options = {}
        if 'spatial' in file_options:
            for k, v in file_options['spatial'].items():
                values_path = os.path.join(file_path, v)
                if check_exist and not os.path.exists(values_path):
                    continue
                options[k] = v

        if 'mappings' in file_options:
            for k, v in file_options['mappings'].items():
                k_name = k.split('/')[-1]
                values_path = os.path.join(file_path, k)
                if check_exist and not os.path.exists(values_path):
                    continue
                m_key = k_name.upper()
                mapping = { "key": k }
                if v is not None and v != 'null':
                    mapping["dims"] = v
                options.setdefault('mappings', {})[m_key] = mapping
        
        if 'factors' in file_options:
            for factor in file_options['factors']:
                values_path = os.path.join(file_path, factor)
                if check_exist and not os.path.exists(values_path):
                    continue    
                options.setdefault('factors', []).append( factor )

    elif file_type == 'anndata-cell-sets.zarr':
        options = []
        if 'sets' in file_options:
            for cell_set in file_options['sets']:
                cell_set_name = cell_set.split('/')[-1]
                if check_exist and not os.path.exists(os.path.join(file_path, cell_set)):
                    continue    
                options.append({
                    "groupName": "".join(w.capitalize() for w in cell_set_name.split("_")),
                    "setName": cell_set
                })

    elif file_type == 'anndata-expression-matrix.zarr':
        options = {}
        for k, v in file_options.items():
            if check_exist and not os.path.exists(os.path.join(file_path, v)):
                continue
            options[k]= v
    
    return options


def write_raster_json(
    title='',
    dataset='',
    file_paths=None,
    files_dir='',
    zarr_dirs=None,
    url='',
    outdir='',
    config_filename='config.json',
    options={}
    ):

    config_json = {
        "name": title,
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
        for file_name, file_type in DATA_TYPES[data_type]:
            
            # first file found will be used in the config file
            file_exists = False
            if file_paths is not None:
                if file_name in [ os.path.basename(x) for x in file_paths ]:
                    idx = [ os.path.basename(x) for x in file_paths ].index(file_name)
                    file_path = file_paths[idx]
                    file_exists = True
            else:
                file_path = os.path.join(files_dir, file_name)
                file_exists = os.path.exists(file_path)

            if file_exists:
                file_entry = {
                    "url": os.path.join(url, file_name),
                    "type": data_type,
                    "fileType": file_type
                }
                file_options = build_options(file_type, file_path, options[file_type] if file_type in options else None)
                if file_options is not None:
                    file_entry["options"] = file_options
                files.append(file_entry)
                break

    # TODO: add image zarr data
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