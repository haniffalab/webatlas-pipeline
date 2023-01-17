#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

// Default params
params.max_n_worker = 30

params.title = ""
params.dataset = ""

params.tsv_delimiter = ","

params.images = []
params.args = []
params.url = ""

params.layout = "minimal"
params.custom_layout = ""
params.outdir = ""
params.config_files = []
params.options = [:]

params.config_params = [
    url: params.url,
    options: params.options,
    layout: params.layout,
    custom_layout: params.custom_layout
]

// if directly writing to s3
params.s3 = false
params.s3_keys = [
    "YOUR_ACCESS_KEY",
    "YOUR_SECRET_KEY"
]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"


verbose_log = true
version = "0.0.1"
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"

process image_to_zarr {
    tag "${image}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-image-to-zarr:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(img_type), path(image)

    output:
    tuple val(stem), val(img_type), path("${stem_str}-${img_type}.zarr"), emit: img_zarr
    tuple val(stem), val(img_type), path("${stem_str}-${img_type}.zarr/OME/METADATA.ome.xml"), emit: ome_xml

    script:
    stem_str = stem.join("-")
    """
    if tiffinfo ${image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        tiffcp -c none ${image} uncompressed.tif
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs uncompressed.tif ${stem_str}-${img_type}.zarr
    else
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs ${image} ${stem_str}-${img_type}.zarr
    fi
    consolidate_md.py ${stem_str}-${img_type}.zarr
    """
}

process ome_zarr_metadata{
    tag "${zarr}"
    debug verbose_log
    container "haniffalab/vitessce-pipeline-processing:${version}"

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

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem_str}*"), emit: converted_files, optional: true

    script:
    stem_str = stem.join("-")
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --path ${file} --stem ${stem_str} ${args_str}
    """
}

process Build_config {
    tag "${stem}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-build-config:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(files), val(img_map), val(config_map)

    output:
    path("${stem_str}-config.json")

    script:
    stem_str = stem.join("-")
    file_paths = files.collect{ /"/ + it + /"/ }.join(",")
    imgs_str = img_map ? "--images '" + new JsonBuilder(img_map).toString() + "'" : ""
    url_str = config_map.url?.trim() ? "--url ${config_map.url.replaceFirst(/\/*$/, "")}/${version}" : ""
    options_str = config_map.options ? "--options '" + (config_map.options instanceof String ? options : new JsonBuilder(config_map.options).toString()) + "'" : ""
    clayout_str = config_map.custom_layout?.trim() ? "--custom_layout \"${config_map.custom_layout}\"" : ""
    """
    build_config.py \
        --title "${stem[0]}" \
        --dataset ${stem[1]} \
        --file_paths '[${file_paths}]' \
        ${imgs_str} \
        ${url_str} \
        ${options_str} \
        --layout ${config_map.layout} ${clayout_str}
    """
}

process Generate_label_image {
    tag "${stem}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), path(file_path), val(file_type), path(ref_img), val(args)

    output:
    tuple val(stem), val("label"), path("${stem_str}-label.tif")

    script:
    stem_str = stem.join("-")
    ref_img_str = ref_img.name != "NO_REF" ? "--ref_img ${ref_img}" : ""
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    generate_label.py --stem ${stem_str} --file_type ${file_type} --file_path ${file_path} ${ref_img_str} ${args_str}
    """
}


Channel.fromPath(params.tsv)
    .splitCsv(header:true, sep:params.tsv_delimiter, quote:"'")
    .map { l -> tuple( tuple(l.title, l.dataset), l ) }
    .branch { stem, l ->
        data: l.data_type in ["h5ad","spaceranger","molecules"]
        images: l.data_type in ["raw_image","label_image","label_image_data"]
        config_params: l.data_type in ["description","url","options"]
            return [stem, ["${l.data_type}": l.data_path]]
        other: true
    }
    .set{inputs}


workflow Full_pipeline {

    Process_files()

    Process_images()


    // Map workflows' outputs to:
    // tuple val(stem), val(files), val(img_map), val(config_map)

    Process_images.out.img_zarrs
        .branch { stem, type, img ->
            raw: type == "raw"
            label: type == "label"
        }
        .set{img_zarrs}

    img_zarrs.raw
        .join(img_zarrs.label)
        .map { stem, raw_type, raw_imgs, label_type, label_imgs -> [
            stem, [raw: raw_imgs, label: label_imgs]
        ]}
        .set{img_map}
    
    inputs.config_params
        .groupTuple()
        .map { stem, it ->
            [ stem, params.config_params << it.collectEntries() { 
                i -> i.collectEntries { k, v -> [(k.toString()): v] } 
                }.findAll { it.value?.trim() ? true : false }
            ]
        }
        .set{config_map}

    Process_files.out.file_paths
        .join(img_map)
        .join(config_map)
        .set{data_for_config}
    

    Build_config(
        data_for_config
        )
}

workflow Process_files {
    if (inputs.data){
        // Map inputs to: 
        // tuple val(stem), path(file), val(type), val(args)
        data_list = inputs.data.flatMap { stem, data_map ->
            data_map.data_path ? 
            [
                [
                    stem,
                    data_map.data_path,
                    data_map.data_type,
                    (data_map.args && data_map.args?.trim() ?
                        data_map.args?.trim() :
                        params.args[data_map.data_type] ?: [:]
                    )
                ]
            ] : [:]
        }

        route_file(data_list)

        files = route_file.out.converted_files.groupTuple(by:0)
        file_paths = files.map { stem, it -> 
            [ stem, it.name ]
        }
    } else {
        files = []
        file_paths = []
    }

    emit:
    files = files
    file_paths = file_paths
}

workflow Process_images {
    if (inputs.images) {
        // Map tif inputs to: 
        // tuple val(stem), val(img_type), path(image)
        img_tifs = inputs.images.filter { stem, data_map ->
            data_map.data_type in ["raw_image", "label_image"]
        }
        .map { stem, data_map ->
            [ 
                stem,
                data_map.data_type.replace("_image",""),
                data_map.data_path
            ]
        }

        // Map label data inputs to: 
        // tuple val(stem), path(file_path), val(file_type), path(ref_img), val(args)
        // Pop file_type (required) and ref_img (optional) from args
        img_data = inputs.images.filter { stem, data_map ->
            data_map.data_type == "label_image_data"
        }
        .map { stem, data_map ->
            [
                stem,
                data_map.data_path,
                *[
                    new JsonSlurper().parseText(data_map.args).file_type,
                    new JsonSlurper().parseText(data_map.args).containsKey("ref_img") ? 
                        new JsonSlurper().parseText(data_map.args).ref_img :
                        file("NO_REF"),
                    new JsonSlurper().parseText(data_map.args)
                        .findAll { !(it.key in ["file_type", "ref_img"]) }
                ]
            ]
        }

        Generate_label_image(img_data)

        all_tifs = img_tifs.mix(Generate_label_image.out)
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
    } else {
        img_zarrs = []
    }

    emit:
    img_zarrs = img_zarrs
}
