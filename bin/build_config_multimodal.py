#!/usr/bin/env python3
"""
build_config.py
====================================
Generates a Vitessce View config
"""

from __future__ import annotations
import typing as T
import os
import fire
import json
from vitessce import (
    VitessceConfig,
    DataType as dt,
    FileType as ft,
    CoordinationType as ct,
    ViewType as vt,
)
from constants import ANNDATA_ZARR_SUFFIX, MOLECULES_JSON_SUFFIX


def write_json(
    project: str = "",
    datasets: dict[str, dict[str]] = {},
    url: str = "",
    integrated: bool = False,
    layout: str = "minimal",
    custom_layout: str = None,
    config_filename_suffix: str = "config.json",
    title: str = "",
    description: str = "",
    outdir: str = "./",
) -> None:
    """This function writes a Vitessce View config JSON file

    Args:
        project (str, optional): Project name. Defaults to "".
        datasets (dict[str, dict[str]], optional): Dictionary of datasets.
            Expected structure: { dataset_name: {
                "file_paths" : [],
                "images": {"raw": [], "label": []},
                "options": {},
                "extended_features": [],
                "obs_type": "cell"
                "is_spatial": True // if has images should be enough
                }
            }
            Defaults to {}.
        url (str, optional): URL to prepend to each file in the config file.
            The URL to the local or remote server that will serve the files.
            Defaults to "".
        integrated (bool, optional): Whether datasets match on features and should
            share same coordination values on them.
        layout (str, optional): Type of predefined layout to use. Defaults to "minimal".
        custom_layout (str, optional): String defining a Vitessce layout following its
            alternative syntax.
            https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout
            https://github.com/vitessce/vitessce-python/blob/1e100e4f3f6b2389a899552dffe90716ffafc6d5/vitessce/config.py#L855
            Defaults to None.
        config_filename_suffix (str, optional): Config filename suffix.
            Defaults to "config.json".
        outdir (str, optional): Directory in which the config file will be written to.
            Defaults to "./".
    """

    config = VitessceConfig(
        "1.0.15",
        name=str(title) if len(title) else str(project),
        description=description,
    )

    for dataset_name in datasets.keys():
        dataset = datasets[dataset_name]
        config_dataset = config.add_dataset(name=dataset_name, uid=dataset_name)
        has_multiple_features = False
        # Add non-image files
        for file_path in dataset["file_paths"]:
            if ANNDATA_ZARR_SUFFIX in file_path:
                # if 1 featureType (e.g. genes) then 1 OBS_FEATURE_MATRIX_ANNDATA_ZARR
                # if 2+ featureTypes (e.g. genes + celltypes) then 1 + n featureTypes OBS_FEATURE_MATRIX_ANNDATA_ZARR
                if "extended_features" in dataset and len(dataset["extended_features"]):
                    has_multiple_features = True

                    if isinstance(dataset["extended_features"], str):
                        dataset["extended_features"] = [dataset["extended_features"]]
                    # Add 'combined' matrix
                    config_dataset.add_file(
                        ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR,
                        url=os.path.join(url, os.path.basename(file_path)),
                        options={"path": dataset["options"]["matrix"] or "X"},
                        coordination_values={
                            ct.OBS_TYPE.value: dataset["obs_type"]
                            if "obs_type" in dataset
                            else "cell",
                            ct.FEATURE_TYPE.value: "combined",
                            ct.FEATURE_VALUE_TYPE.value: "expression",
                        },
                    )

                    # Add gene matrix
                    config_dataset.add_file(
                        ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR,
                        url=os.path.join(url, os.path.basename(file_path)),
                        options={
                            "path": dataset["options"]["matrix"] or "X",
                            "featureFilterPath": "var/is_gene",
                        },
                        coordination_values={
                            ct.OBS_TYPE.value: dataset["obs_type"]
                            if "obs_type" in dataset
                            else "cell",
                            ct.FEATURE_TYPE.value: "gene",
                            ct.FEATURE_VALUE_TYPE.value: "expression",
                        },
                    )

                    for feature in dataset["extended_features"]:
                        # Add extended matrix
                        config_dataset.add_file(
                            ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR,
                            url=os.path.join(url, os.path.basename(file_path)),
                            options={
                                "path": dataset["options"]["matrix"] or "X",
                                "featureFilterPath": f"var/is_{feature}",
                            },
                            coordination_values={
                                ct.OBS_TYPE.value: dataset["obs_type"]
                                if "obs_type" in dataset
                                else "cell",
                                ct.FEATURE_TYPE.value: dataset,
                                ct.FEATURE_TYPE.value: "abundance"
                                if feature == "celltype"
                                else "expression",
                            },
                        )

                    # Add other options (locations, embeddings, labels)
                    options = build_anndatazarr_options(
                        dataset["options"], include_matrix=False
                    )
                    if len(options.keys()):
                        config_dataset.add_file(
                            ft.ANNDATA_ZARR,
                            url=os.path.join(url, os.path.basename(file_path)),
                            options=options,
                            coordination_values={
                                ct.OBS_TYPE.value: dataset["obs_type"]
                                if "obs_type" in dataset
                                else "cell"
                            },
                        )

                else:
                    config_dataset.add_file(
                        ft.ANNDATA_ZARR,
                        url=os.path.join(url, os.path.basename(file_path)),
                        options=build_anndatazarr_options(dataset["options"]),
                        coordination_values={
                            ct.OBS_TYPE.value: dataset["obs_type"]
                            if "obs_type" in dataset
                            else "cell"
                        },
                    )

            if MOLECULES_JSON_SUFFIX in file_path:
                config_dataset.add_file(
                    ft.MOLECULES_JSON,
                    url=os.path.join(url, os.path.basename(file_path)),
                )

        # Add image files if any
        # if (
        #     "images" in dataset
        #     and dataset["images"].keys()
        #     and any([len(dataset["images"][k]) for k in dataset["images"].keys()])
        # ):
        if dataset["is_spatial"]:
            config_dataset.add_file(
                ft.RASTER_JSON,
                options=build_raster_options(dataset["images"], url=url),
                coordination_values={
                    ct.OBS_TYPE.value: dataset["obs_type"]
                    if "obs_type" in dataset
                    else "cell",
                    ct.FEATURE_TYPE.value: "combined"
                    if has_multiple_features
                    else "expression",
                },
            )

        config_json = config.to_dict()

        if outdir and not os.path.isdir(outdir):
            os.mkdir(outdir)
        with open(
            os.path.join(
                outdir or "", f"{project}-multimodal-{config_filename_suffix}"
            ),
            "w",
        ) as out_file:
            json.dump(config_json, out_file, indent=2)


def build_anndatazarr_options(
    file_options: dict[str, T.Any], include_matrix: bool = True
):
    # @TODO: Change yaml input keys to 'locations', 'embeddings', 'labels'
    options = {}
    if "spatial" in file_options and "xy" in file_options["spatial"]:
        options[dt.OBS_LOCATIONS.value] = {"path": file_options["spatial"]["xy"]}
    if "mappings" in file_options:
        for k, v in file_options["mappings"].items():
            embedding_type = k.split("/")[-1].upper()
            v = v if v is not None and v != "null" else [0, 1]
            embedding = {"path": k, "embeddingType": embedding_type, "dims": v}
            options.setdefault(dt.OBS_EMBEDDING.value, []).append(embedding)
    if "factors" in file_options:
        for factor in file_options["factors"]:
            label_type = factor.split("/")[-1].capitalize()
            label = {"path": factor, "obsLabelsType": label_type}
            options.setdefault(dt.OBS_LABELS.value, []).append(label)
    if "sets" in file_options:
        for obs_set in file_options["sets"]:
            if type(obs_set) != dict:
                obs_set = {"name": obs_set}
            obs_set_options = {
                "name": "".join(
                    w.capitalize() for w in obs_set["name"].split("/")[-1].split("_")
                ),
                "path": obs_set["name"],
            }
            if "score" in obs_set:
                obs_set_options["name"] += " with scores"
                obs_set_options["scorePath"] = obs_set["score"]
            options.setdefault(dt.OBS_SETS.value, []).append(obs_set_options)
    if include_matrix and "matrix" in file_options:
        options[dt.OBS_FEATURE_MATRIX.value] = {"path": file_options["matrix"]}

    return options


def build_raster_options(
    images: dict[str, list[dict[str, T.Any]]], url: str
) -> dict[str, T.Any]:
    """Function that creates the View config's options for image files

    Args:
        images (dict[str, list[dict[str, T.Any]]], optional): Dictionary containing for each image type key (raw and label)
            a list of dictionaries (one per image of that type) with the corresponding path and metadata for that image.
            Defaults to {}.
        url (str): URL to prepend to each file in the config file.
            The URL to the local or remote server that will serve the files

    Returns:
        dict[str, T.Any]: Options dictionary for View config file
    """
    raster_options = {"renderLayers": [], "schemaVersion": "0.0.2", "images": []}
    for img_type in images.keys():  # "raw" or "label"
        for img in images[img_type]:
            image_name = os.path.splitext(os.path.basename(img["path"]))[0]
            channel_names = (
                img["md"]["channel_names"]
                if "channel_names" in img["md"] and len(img["md"]["channel_names"])
                else (
                    ["Labels"]
                    if img_type == "label"
                    else [f"Channel {x}" for x in range(int(img["md"]["C"]))]
                )
            )

            raster_options["renderLayers"].append(image_name)
            raster_options["images"].append(
                {
                    "name": image_name,
                    "url": os.path.join(url, os.path.basename(img["path"])),
                    "type": "zarr",
                    "metadata": {
                        "isBitmask": img_type == "label",
                        "dimensions": [
                            {"field": "t", "type": "quantitative", "values": None},
                            {
                                "field": "channel",
                                "type": "nominal",
                                "values": channel_names,
                            },
                            {"field": "y", "type": "quantitative", "values": None},
                            {"field": "x", "type": "quantitative", "values": None},
                        ],
                        "isPyramid": True,
                        "transform": {"translate": {"y": 0, "x": 0}, "scale": 1},
                    },
                }
            )

    return raster_options


if __name__ == "__main__":
    fire.Fire(write_json)
