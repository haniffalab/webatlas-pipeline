#!/usr/bin/env python3
"""
build_config.py
====================================
Generates a Vitessce View config
"""

from __future__ import annotations
from itertools import zip_longest
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
    hconcat,
    vconcat,
)
from constants.suffixes import ANNDATA_ZARR_SUFFIX, MOLECULES_JSON_SUFFIX


def write_json(
    project: str = "",
    datasets: dict[str, dict[str]] = {},
    extended_features: T.Union[str, list] = [],
    url: str = "",
    config_filename_suffix: str = "config.json",
    title: str = "",
    description: str = "",
    outdir: str = "./",
) -> None:
    """This function writes a Vitessce View config JSON file

    Args:
        project (str, optional): Project name. Defaults to "".
        datasets (dict[str, dict[str]], optional): Dictionary of datasets.
            Expected structure: { dataset_name: { "file_paths" : [],
            "images": {"raw": [], "label": []},
            "options": {}, "obs_type": "cell",
            "is_spatial": True } }
            Defaults to {}.
        extended_features (Union[list[str], str], optional): List of features or
            string of single feature on which the expression matrix was extended
            and var/is_{feature} is present.
            Defaults to [].
        url (str, optional): URL to prepend to each file in the config file.
            The URL to the local or remote server that will serve the files.
            Defaults to "".
        config_filename_suffix (str, optional): Config filename suffix.
            Defaults to "config.json".
        title (str, optional): Data title to show in the visualization. Defaults to "".
        description (str, optional): Data description to show in the visualization.
            Defaults to "".
        outdir (str, optional): Directory in which the config file will be written to.
            Defaults to "./".
    """
    config = VitessceConfig(
        "1.0.15",
        name=str(title) if len(title) else str(project),
        description=description,
    )

    features = {
        "gene": {
            "value_type": "expression",
            "ftype": config.add_coordination(ct.FEATURE_TYPE)[0].set_value("gene"),
            "fvalue_type": config.add_coordination(ct.FEATURE_VALUE_TYPE)[0].set_value(
                "expression"
            ),
        }
    }

    has_multiple_features = False
    if len(extended_features):
        has_multiple_features = True
        if isinstance(extended_features, str):
            extended_features = [extended_features]

        # Set coordination types for features
        for feature in extended_features:
            feature_value_type = "abundance" if feature == "celltype" else "expression"
            features.setdefault(feature, {})["value_type"] = feature_value_type
            features.setdefault(feature, {})["ftype"] = config.add_coordination(
                ct.FEATURE_TYPE
            )[0].set_value(feature)
            features.setdefault(feature, {})["fvalue_type"] = config.add_coordination(
                ct.FEATURE_VALUE_TYPE
            )[0].set_value(feature_value_type)

        multi_ftype, *_ = config.add_coordination(ct.FEATURE_TYPE)
        multi_ftype.set_value("combined")

    # Use dicts to save coordination types for observations and embeddings
    obs_coordination = {}
    embeddings_coordination = {}

    config_views = []

    # Add datasets to config
    for dataset_name in datasets.keys():
        dataset = datasets[dataset_name]
        config_dataset = config.add_dataset(name=dataset_name, uid=dataset_name)

        if isinstance(dataset["is_spatial"], str):
            dataset["is_spatial"] = dataset["is_spatial"].lower() == "true"
        dataset["is_spatial"] = (
            dataset["is_spatial"] and dataset.get("images", {}).keys()
        )

        dataset_obs_type = dataset.get("obs_type", "cell")

        if dataset_obs_type not in obs_coordination:
            obs_coordination[dataset_obs_type] = config.add_coordination(ct.OBS_TYPE)[
                0
            ].set_value(dataset_obs_type)

        # Add non-image files
        for file_path in dataset["file_paths"]:
            if ANNDATA_ZARR_SUFFIX in file_path:
                # if 1 featureType (e.g. genes) then 1 OBS_FEATURE_MATRIX_ANNDATA_ZARR
                # if 2+ featureTypes (e.g. genes + celltypes) then 1 + n * OBS_FEATURE_MATRIX_ANNDATA_ZARR
                if has_multiple_features:
                    # Add 'combined' matrix
                    config_dataset.add_file(
                        ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR,
                        url=os.path.join(url, os.path.basename(file_path)),
                        options={"path": dataset.get("options", {}).get("matrix", "X")},
                        coordination_values={
                            ct.OBS_TYPE.value: dataset_obs_type,
                            ct.FEATURE_TYPE.value: "combined",
                            ct.FEATURE_VALUE_TYPE.value: "expression",
                        },
                    )

                    for feature in ["gene", *extended_features]:
                        # Add extended matrix
                        config_dataset.add_file(
                            ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR,
                            url=os.path.join(url, os.path.basename(file_path)),
                            options={
                                "path": dataset.get("options", {}).get("matrix", "X"),
                                "featureFilterPath": f"var/is_{feature}",
                            },
                            coordination_values={
                                ct.OBS_TYPE.value: dataset_obs_type,
                                ct.FEATURE_TYPE.value: feature,
                                ct.FEATURE_VALUE_TYPE.value: features[feature][
                                    "value_type"
                                ],
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
                            coordination_values={ct.OBS_TYPE.value: dataset_obs_type},
                        )

                else:  # just one feature type
                    config_dataset.add_file(
                        ft.ANNDATA_ZARR,
                        url=os.path.join(url, os.path.basename(file_path)),
                        options=build_anndatazarr_options(dataset["options"]),
                        coordination_values={
                            ct.OBS_TYPE.value: dataset_obs_type,
                            ct.FEATURE_TYPE.value: "gene",
                            ct.FEATURE_VALUE_TYPE.value: "expression",
                        },
                    )

                # If not spatial dataset add embedding view
                if not dataset["is_spatial"] and len(
                    dataset.get("options", {}).get("mappings", {})
                ):
                    dataset_embedding_type = (
                        list(dataset["options"]["mappings"].keys())[0]
                        .split("/")[-1]
                        .upper()
                    )
                    if dataset_embedding_type not in embeddings_coordination:
                        (
                            embeddings_coordination[dataset_embedding_type]
                        ) = config.add_coordination(ct.EMBEDDING_TYPE)[0].set_value(
                            dataset_embedding_type
                        )

                    config_views.append(
                        config.add_view(
                            vt.SCATTERPLOT, dataset=config_dataset
                        ).use_coordination(
                            embeddings_coordination[dataset_embedding_type],
                            obs_coordination[dataset_obs_type],
                            multi_ftype
                            if has_multiple_features
                            else features["gene"]["ftype"],
                        )
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
                    ct.OBS_TYPE.value: dataset_obs_type,
                    ct.FEATURE_TYPE.value: "combined"
                    if has_multiple_features
                    else "expression",
                },
            )

            # If spatial dataset add spatial view (no embedding)
            spatial_coord = [
                obs_coordination[dataset_obs_type],
                *config.add_coordination(
                    ct.SPATIAL_IMAGE_LAYER,
                    ct.SPATIAL_SEGMENTATION_LAYER,
                    ct.SPATIAL_ZOOM,
                    ct.SPATIAL_TARGET_X,
                    ct.SPATIAL_TARGET_Y,
                ),
            ]

            config_views.append(
                config.add_view(vt.SPATIAL, dataset=config_dataset).use_coordination(
                    *spatial_coord,
                    multi_ftype if has_multiple_features else features["gene"]["ftype"],
                )
            )
            # Make layer controller optional?
            config_views.append(
                config.add_view(
                    vt.LAYER_CONTROLLER,
                    dataset=config_dataset,
                ).use_coordination(*spatial_coord)
            )

    # Add feature lists
    for feature in features:
        config_views.append(
            config.add_view(
                vt.FEATURE_LIST,
                dataset_uid=list(datasets.keys())[0],
            ).use_coordination(
                features[feature]["ftype"], features[feature]["fvalue_type"]
            )
        )

    # Set layout
    config.layout(concat_views(config_views))

    # Write json
    config_json = config.to_dict()
    # Ensure views' grid coorindates are integers
    for component in config_json["layout"]:
        for k in ("x", "y", "w", "h"):
            component[k] = round(component[k])

    if outdir and not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(
        os.path.join(outdir or "", f"{project}-multimodal-{config_filename_suffix}"),
        "w",
    ) as out_file:
        json.dump(config_json, out_file, indent=2)


def concat_views(views: list, axis: str = "v"):
    """Recursively concatenate views"""
    if len(views) <= 2:
        return hconcat(*views) if axis == "h" else vconcat(*views)
    nviews = []
    # Make pairs and concatenate them
    for i, j in zip_longest(views[::2], views[1::2]):
        if j is not None:
            nviews.append(hconcat(i, j) if axis == "h" else vconcat(i, j))
        else:
            nviews.append(i)
    # Recursive to concatenate pairs of concatenated pairs
    return concat_views(nviews, "v" if axis == "h" else "h")


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
        images (dict[str, list[dict[str, T.Any]]], optional): Dictionary containing
            for each image type key (raw and label) a list of dictionaries
            (one per image of that type) with the corresponding path and
            metadata for that image.
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
                if len(img.get("md", {}).get("channel_names", []))
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
