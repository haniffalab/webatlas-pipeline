#!/usr/bin/env python3

from collections import defaultdict
import os
import fire
import json
import re
from itertools import chain, cycle

from xml.etree import ElementTree as ET
import pandas as pd

from vitessce import (
    VitessceConfig,
    FileType as ft,
    CoordinationType as ct,
    Component as cm
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
    layout='minimal',
    custom_layout=None,
    codebook='',
    ):

    print(f"zarr paths : {zarr_dirs}")

    NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

    ome_metadata = ET.parse("raw_image/OME/METADATA.ome.xml")
    dimOrder = ome_metadata.find("./*/ome:Pixels", NS).attrib["DimensionOrder"]
    print(dimOrder)
    X = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeX"]
    Y = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeY"]
    Z = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeZ"]
    C = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeC"]
    T = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeT"]
    print(X, Y, Z, C, T)
    channels = [channel.attrib["Name"] for channel in ome_metadata.findall("./**/ome:Channel", NS) if "Name" in channel.attrib]
    print(channels)

    codebook = pd.read_csv(codebook)
    print(codebook)

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
                # Set a coordination scope for any 'mapping'
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
    # Set layout with alternative syntax https://github.com/vitessce/vitessce-python/blob/1e100e4f3f6b2389a899552dffe90716ffafc6d5/vitessce/config.py#L855
    config_layout = custom_layout if custom_layout and len(custom_layout) else DEFAULT_LAYOUTS[layout]
    views, views_indx = [], []
    for m in re.finditer("[a-zA-Z]+", config_layout):
        # TODO: catch error
        component = cm(m.group())

        # Remove component from layout if its required data type is not present
        if component in COMPONENTS_DATA_TYPES and COMPONENTS_DATA_TYPES[component].isdisjoint(dts):
            # Remove trailing or leading operators along with component
            trail = m.end() < len(config_layout) and config_layout[m.end()] in ["|","/"]
            lead = not trail and m.start() > 0 and config_layout[m.start()-1] in ["|","/"]
            views_indx.append((m.start()-1*lead, m.end()+1*trail, ""))
            continue

        view = config.add_view(component, dataset=config_dataset)

        if component in COMPONENTS_COORDINATION_TYPES:
            for coordination_type in COMPONENTS_COORDINATION_TYPES[component]:
                # Cycle through coordination scopes with the required coordination type
                # TODO: catch error
                view.use_coordination(next(coordination_types[coordination_type]))

        views.append(view)
        views_indx.append((m.start(), m.end(), "views[{}]".format(len(views)-1)))

    # Replace components names in layout string with the corresponding view object
    for (start, end, view) in sorted(views_indx, key=lambda x: x[0], reverse=True):
        config_layout = config_layout[:start] + view + config_layout[end:]

    # Concatenate views
    exec("config.layout({})".format(config_layout))

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
