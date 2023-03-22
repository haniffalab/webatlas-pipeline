#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

verbose_log = true
version = "0.0.1"

//////////////////////////////////////////////////////

// Default params
params.max_n_worker = 30

params.outdir = ""

params.vitessce_options = [:]
params.layout = "minimal"
params.custom_layout = ""

params.vitessce_config_map = [
    url: "http://localhost:3000/",
    options: params.vitessce_options,
    layout: params.layout,
    custom_layout: params.custom_layout,
    title: "",
    description: ""
]

// if directly writing to s3
params.s3 = false
params.s3_keys = [
    "YOUR_ACCESS_KEY",
    "YOUR_SECRET_KEY"
]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"

//////////////////////////////////////////////////////

data_types = ["h5ad","spaceranger","xenium","merscope","molecules"]
image_types = ["raw_image","label_image","raw_image_data","label_image_data"]
vitessce_params = ["title","description","url","vitessce_options","layout","custom_layout"]

outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"

//////////////////////////////////////////////////////

// Get inputs
Channel.from(params.projects)
    .map { p -> [p.project, p.datasets] }
    .transpose()
    .multiMap {
        project, dataset -> 
            data: [ tuple(project, dataset.dataset), dataset.data ]
            config_map: [
                tuple(project, dataset.dataset),
                params.vitessce_config_map + dataset.subMap(vitessce_params)
            ]
    }
    .set {datasets}

datasets.data
    .transpose(by:1)
    .branch{ stem, d ->
        data: d.data_type in data_types
        images: d.data_type in image_types
        other: true
    }
    .set{inputs}

inputs.other
    .collect { stem, d -> d.data_type }
    .view{ "Unrecognized data_type(s) ${it.unique()}" }

//////////////////

interm_dt = [
    spaceranger: ["h5ad"],
    xenium: ["h5ad"],
    merscope: ["h5ad"]
]

project_args = [:]
dataset_args = [:]

params.projects.each{ p ->
    project_args[p.project] = p.args ?: [:]
    p.datasets.each{ d ->
        dataset_args[[p.project, d.dataset]] = d.args ?: [:]
    }
}

def getSubMapValues (m, keys) {
    m.subMap(keys).values().sum() ?: [:]
}

def mergeArgs (stem, data_type, args) {
    getSubMapValues(params.args, [data_type, *interm_dt[data_type]]) + 
    getSubMapValues(project_args[stem[0]], [data_type, *interm_dt[data_type]]) + 
    getSubMapValues(dataset_args[stem], [data_type, *interm_dt[data_type]]) + 
    (args ?: [:])
}

//////////////////////////////////////////////////////

process image_to_zarr {
    tag "${image}"
    debug verbose_log

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(prefix), val(img_type), path(image), val(keep_filename)

    output:
    tuple val(stem), val(img_type), path("${filename}.zarr"), emit: img_zarr
    tuple val(stem), val(img_type), path("${filename}.zarr/OME/METADATA.ome.xml"), emit: ome_xml

    script:
    filename = keep_filename ? image.baseName : ([*stem, prefix, img_type] - null - "").join("-")
    """
    if tiffinfo ${image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        tiffcp -c none ${image} uncompressed.tif
        bioformats2raw --no-hcs uncompressed.tif ${filename}.zarr
    else
        bioformats2raw --no-hcs ${image} ${filename}.zarr
    fi
    consolidate_md.py ${filename}.zarr
    """
}

process ome_zarr_metadata{
    tag "${zarr}"
    debug verbose_log

    input:
    tuple val(stem), val(img_type), path(zarr)

    output:
    tuple val(stem), val(img_type), stdout

    script:
    """
    ome_zarr_metadata.py --xml_path ${zarr}
    """
}

process route_file {
    tag "${type}, ${file}"
    debug verbose_log
    cache "lenient"

    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), val(prefix), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem_str}*"), emit: converted_files, optional: true
    tuple val(stem), path("tmp-${stem_str}*"), emit: extra_files, optional: true

    script:
    stem_str = ([*stem, prefix] - null - "").join("-")
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --path ${file} --stem ${stem_str} ${args_str}
    """
}

process Build_config {
    tag "${stem}"
    debug verbose_log
    cache false

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(config_map), val(files), val(img_map), 

    output:
    path("${stem_str}-config.json")

    script:
    stem_str = stem.join("-")
    file_paths = files.collect{ /"/ + it + /"/ }.join(",")
    imgs_str = img_map ? "--images '" + new JsonBuilder(img_map).toString() + "'" : ""
    url_str = config_map.url?.trim() ? "--url \"${config_map.url.trim()}\"" : ""
    options_str = config_map.options ? "--options '" + (config_map.options instanceof String ? options : new JsonBuilder(config_map.options).toString()) + "'" : ""
    clayout_str = config_map.custom_layout?.trim() ? "--custom_layout \"${config_map.custom_layout}\"" : ""
    """
    build_config.py \
        --project "${stem[0]}" \
        --dataset "${stem[1]}" \
        --file_paths '[${file_paths}]' \
        ${imgs_str} \
        ${url_str} \
        ${options_str} \
        --layout "${config_map.layout}" ${clayout_str} \
        --title "${config_map.title}" \
        --description "${config_map.description}"
    """
}

process Generate_image {
    tag "${stem}, ${img_type}, ${file_path}"
    debug verbose_log

    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), val(prefix), val(img_type), path(file_path), val(file_type), path(ref_img), val(args)

    output:
    tuple val(stem), val(prefix), val(img_type), path("${stem_str}*.tif")

    script:
    stem_str = ([*stem, prefix] - null - "").join("-")
    ref_img_str = ref_img.name != "NO_REF" ? "--ref_img ${ref_img}" : ""
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    generate_image.py \
        --stem ${stem_str} \
        --img_type ${img_type} \
        --file_type ${file_type} \
        --file_path ${file_path} \
        ${ref_img_str} ${args_str}
    """
}

//////////////////////////////////////////////////////

workflow Full_pipeline {

    Process_files()

    Process_images()

    Output_to_config(
        Process_files.out.file_paths,
        Process_images.out.img_zarrs
        )
    
}


workflow Process_files {
    // Map inputs to: 
    // tuple val(stem), val(prefix), path(file), val(type), val(args)
    data_list = inputs.data.flatMap { stem, data_map ->
        data_map.data_path ?
        [
            [
                stem,
                data_map.prefix ?: "",
                data_map.data_path,
                data_map.data_type,
                mergeArgs(stem, data_map.data_type, data_map.args)
            ]
        ] : [:]
    }

    route_file(data_list)
    files = route_file.out.converted_files
        .map { stem, paths ->
            [ stem, [paths].flatten() ]
        }
        .transpose(by: 1)
        .groupTuple(by:0)
    file_paths = files.map { stem, it -> 
        [ stem, it.name ]
    }

    emit:
    files = files
    file_paths = file_paths
}


workflow Process_images {
    // Map tif inputs to:
    // tuple val(stem), val(prefix), val(img_type), path(image)
    img_tifs = inputs.images.filter { stem, data_map ->
        data_map.data_type in ["raw_image", "label_image"]
    }
    .map { stem, data_map ->
        [ 
            stem,
            data_map.prefix,
            data_map.data_type.replace("_image",""),
            data_map.data_path,
            false // keep_filename
        ]
    }

    // Map raw/label data inputs to:
    // tuple val(stem), val(prefix), val(img_type), path(file_path), val(file_type), path(ref_img), val(args)
    img_data = inputs.images.filter { stem, data_map ->
        data_map.data_type in ["raw_image_data", "label_image_data"]
    }
    .map { stem, data_map ->
        [
            stem,
            data_map.prefix,
            data_map.data_type.replace("_image_data",""),
            data_map.data_path,
            data_map.file_type,
            data_map.ref_img ?: file("NO_REF"),
            data_map.args ?: [:]
        ]
    }

    Generate_image(img_data)

    Generate_image.out
        .map { stem, prefix, type, paths ->
            [
                stem,
                prefix,
                type,
                [paths].flatten(),
                true // keep_filename
            ]
        }
        .transpose(by: 3)
        .set {label_tifs}

    all_tifs = img_tifs.mix(label_tifs)
    image_to_zarr(all_tifs)

    ome_zarr_metadata(image_to_zarr.out.ome_xml)

    img_zarrs = image_to_zarr.out.img_zarr
        .join(ome_zarr_metadata.out, by: [0,1])
        .map { stem, type, path, md ->
            [
                stem, type, [path: path.name, md: new JsonSlurper().parseText(md.trim())]
            ]
        }
        .groupTuple(by: [0,1])

    emit:
    img_zarrs = img_zarrs
}


workflow Output_to_config {
    take: out_file_paths
    take: out_img_zarrs
    main:

        // Map workflows' outputs to:
        // tuple val(stem), val(files), val(img_map), val(config_map)

        out_img_zarrs
            .map { stem, type, img -> 
                [stem, [type: type, img: img]]
            }
            .branch { stem, data ->
                raw: data.type == "raw"
                label: data.type == "label"
            }
            .set{img_zarrs}

        img_zarrs.raw
            .join(img_zarrs.label, remainder: true)
            .map { stem, raw_data, label_data -> [
                stem,
                [
                    raw: raw_data ? raw_data.img : [],
                    label: label_data ? label_data.img : []
                ]
            ]}
            .set{img_map}

        datasets.config_map
            .join(out_file_paths, remainder: true)
            .join(img_map, remainder: true)
            .set{data_for_config}

        Build_config(
            data_for_config
            )
}