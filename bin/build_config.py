#!/usr/bin/env python3

from collections import defaultdict
import os
import fire
import json
from itertools import chain, cycle

from vitessce import (
    VitessceConfig,
    FileType as ft,
    CoordinationType as ct,
    hconcat, vconcat
)
from constants import (
    DATA_TYPES, 
    DEFAULT_OPTIONS, DEFAULT_LAYOUTS, 
    COMPONENTS_COORDINATION_TYPES, COMPONENTS_DATA_TYPES
)

def build_options(file_type, file_path, file_options=None, check_exist=False):
    options = None
    if file_options is None:
        file_options = DEFAULT_OPTIONS
    
    if file_type == ft.ANNDATA_CELLS_ZARR:
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
                v = v if v is not None and v != 'null' else [0,1]
                mapping["dims"] = v
                options.setdefault('mappings', {})[m_key] = mapping
        
        if 'factors' in file_options:
            for factor in file_options['factors']:
                values_path = os.path.join(file_path, factor)
                if check_exist and not os.path.exists(values_path):
                    continue    
                options.setdefault('factors', []).append( factor )

    elif file_type == ft.ANNDATA_CELL_SETS_ZARR:
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

    elif file_type == ft.ANNDATA_EXPRESSION_MATRIX_ZARR:
        options = {}
        ematrix_ops = set(["matrix","matrixGeneFilter","geneAlias"])
        for k, v in file_options.items():
            if k not in ematrix_ops or (check_exist and not os.path.exists(os.path.join(file_path, v))):
                continue
            options[k]= v
    
    return options


def write_json(
    title='',
    dataset='',
    file_paths=None,
    files_dir='',
    zarr_dirs=None,
    url='',
    outdir='',
    config_filename='config.json',
    options={},
    layout='minimal'
    ):

    config = VitessceConfig()
    config_dataset = config.add_dataset(title, dataset)

    coordination_types = defaultdict(list)
    file_paths_names = { os.path.basename(x):x for x in file_paths or [] }
    dts = set([])
    for data_type in DATA_TYPES:
        for file_name, file_type in DATA_TYPES[data_type]:
            
            # first file type found will be used in the config file
            file_exists = False
            if file_paths is not None:
                if file_name in file_paths_names:
                    file_path = file_paths_names[file_name]
                    file_exists = True
            else:
                file_path = os.path.join(files_dir, file_name)
                file_exists = os.path.exists(file_path)

            if file_exists:
                file_options = build_options(file_type, file_path, options)
                if 'mappings' in file_options:
                    coordination_types[ct.EMBEDDING_TYPE] = cycle(chain(
                        coordination_types[ct.EMBEDDING_TYPE],
                        [config.set_coordination_value(ct.EMBEDDING_TYPE.value, k, k) for k in file_options['mappings']]
                    ))
                config_dataset.add_file(
                    data_type, file_type,
                    url = os.path.join(url, file_name),
                    options = file_options
                )
                dts.add(data_type)
                break
    
    
    # Get layout components/views
    config_layout = []
    for col in DEFAULT_LAYOUTS[layout]:
        new_col = []
        for component in col:
            if component in COMPONENTS_DATA_TYPES and COMPONENTS_DATA_TYPES[component].isdisjoint(dts):
                continue
            view = config.add_view(component, dataset=config_dataset)

            if component in COMPONENTS_COORDINATION_TYPES:
                for coordination_type in COMPONENTS_COORDINATION_TYPES[component]:
                    # TODO: catch error
                    view.use_coordination(next(coordination_types[coordination_type]))

            new_col.append(view)
        if len(new_col): config_layout.append(new_col)
    
    # Concatenate views
    config.layout(hconcat(*[vconcat(*[view for view in col]) for col in config_layout]))
    
    # Make sure views' grid coordinates are integers
    config_json = config.to_dict()
    for l in config_json["layout"]:
        for k in ("x","y","w","h"):
            l[k] = round(l[k])
            

    if outdir and not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir or '', config_filename), "w") as out_file:
        json.dump(config_json, out_file, indent=2)


if __name__ == "__main__":
    fire.Fire(write_json)