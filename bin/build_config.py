#!/usr/bin/env python3

from collections import defaultdict
import os
import fire
import json
import logging
from itertools import chain, cycle

from xml.etree import ElementTree as ET
import pandas as pd


def get_image_basic_metadata(xml_path):
    NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

    ome_metadata = ET.parse(xml_path)
    dimOrder = ome_metadata.find("./*/ome:Pixels", NS).attrib["DimensionOrder"]
    X = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeX"]
    Y = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeY"]
    Z = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeZ"]
    C = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeC"]
    T = ome_metadata.find("./*/ome:Pixels", NS).attrib["SizeT"]
    channels = [channel.attrib["Name"] for channel in ome_metadata.findall("./**/ome:Channel", NS) if "Name" in channel.attrib]
    return dimOrder, channels, X, Y, Z, C, T


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


def build_raster_options(image_zarr, url):
    raster_options = {
        "renderLayers": [],
        "schemaVersion": "0.0.2",
        "images": []
    }
    for image in image_zarr.keys():
            raster_options["renderLayers"].append(image)
            raster_options["images"].append({
                "name": image,
                "url": os.path.join(url, image),
                "type": "zarr",
                "metadata": {
                    "isBitmask": True,
                    "dimensions": [
                        {"field": "t",
                          "type": "quantitative",
                          "values": None },
                        {"field": "channel",
                          "type": "nominal",
                          "values": image_zarr[image]["channel_names"] },
                        {"field": "y",
                          "type": "quantitative",
                          "values": None },
                        {"field": "x",
                          "type": "quantitative",
                          "values": None }
                    ],
                    "isPyramid": True,
                        "transform": {
                            "translate": { "y": 0, "x": 0 },
                        "scale": 1
                    }
                }
            })
    return raster_options


def write_json(
    title='',
    dataset='',
    file_paths=[],
    image_zarr={},
    url='',
    outdir='',
    config_filename='config.json',
    options={},
    layout='minimal',
    custom_layout=None
    ):
    import regex
    from vitessce import (
        VitessceConfig,
        DataType as dt,
        FileType as ft,
        CoordinationType as ct,
        Component as cm
    )
    from constants import (
        DATA_TYPES,
        DEFAULT_OPTIONS, DEFAULT_LAYOUTS,
        COMPONENTS_COORDINATION_TYPES, COMPONENTS_DATA_TYPES
    )

    has_files = False

    config = VitessceConfig()
    config_dataset = config.add_dataset(title, dataset)

    coordination_types = defaultdict(list)
    file_paths_names = { os.path.basename(x):x for x in file_paths }
    dts = set([])

    if len(image_zarr.items()):
        has_files = True
        config_dataset.add_file(
            dt.RASTER, ft.RASTER_JSON,
            options = build_raster_options(image_zarr, url)
        )
        dts.add(dt.RASTER)

    for data_type in DATA_TYPES:
        for file_name, file_type in DATA_TYPES[data_type]:

            # first file type found will be used in the config file
            file_exists = False
            if file_name in file_paths_names:
                file_path = file_paths_names[file_name]
                file_exists = True

            if file_exists:
                has_files = True
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

    if not has_files:
        logging.error("No files to add to config file")
        quit(1)

    # Get layout components/views
    # Set layout with alternative syntax https://github.com/vitessce/vitessce-python/blob/1e100e4f3f6b2389a899552dffe90716ffafc6d5/vitessce/config.py#L855
    config_layout = custom_layout if custom_layout and len(custom_layout) else DEFAULT_LAYOUTS[layout]
    views, views_indx = [], []
    for m in regex.finditer("[a-zA-Z]+", config_layout):
        # TODO: catch error
        component = cm(m.group())

        # Remove component from layout if its required data type is not present
        if component in COMPONENTS_DATA_TYPES and COMPONENTS_DATA_TYPES[component].isdisjoint(dts):
            views_indx.append((m.start(), m.end(), ""))
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

    # Remove remaining | and / operators and empty parenthesis
    empty_par = r'\((?>[^\(\)a-zA-Z]|(?R))*\)'
    unbalanced_op = r'(?<=\()([\/\|]+)|([\/\|]+)(?=\))|(?<=^)([\/\|]+)|([\/\|]+)(?=$)'
    repeated_op = r'(?<=[\/\|])([\/\|]+)'
    config_layout = regex.sub(empty_par, '', config_layout)
    config_layout = regex.sub(unbalanced_op, '', config_layout)
    config_layout = regex.sub(repeated_op, '', config_layout)

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
