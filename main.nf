#!/usr/bin/env/ nextflow

import groovy.json.*

nextflow.enable.dsl=2

verbose_log = true
version = "0.5.2"

//////////////////////////////////////////////////////

// Default params
params.max_n_worker = 30

params.outdir = ""
params.args = [:]
params.projects = []
params.write_spatialdata = false
params.publish_generated_img = false

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
    (dataset_args[stem]?.rotate ?: [:]) + // rotate: {rotate_degrees: 90, spatial_shape: [width, height]}
    (args ?: [:])
}

//////////////////////////////////////////////////////

def warnParams () {
    if (!workflow.commandLine.contains("-params-file")){
        log.warn "No -params-file provided"
    }
}

//////////////////////////////////////////////////////

// @TODO: move rotation to separate process and run with all_tifs
process image_to_zarr {
    tag "${image}"
    debug verbose_log

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(prefix), val(img_type), path(image), val(keep_filename), val(rotate_degrees)

    output:
    tuple val(stem), val(img_type), path("${filename}.zarr"), emit: img_zarr
    tuple val(stem), val(img_type), path("${filename}.zarr/OME/METADATA.ome.xml"), emit: ome_xml

    script:
    filename = keep_filename ? image.baseName : ([*stem, prefix, img_type] - null - "").join("-")
    tmp_image = "tmp-${filename}.tif"
    """
    if [[ ${rotate_degrees} != "NO_ROT" ]]
    then
        if [[ ${rotate_degrees} != 90 && ${rotate_degrees} != 180 && ${rotate_degrees} != 270 ]]
        then
            echo "Invalid rotation value: ${rotate_degrees}"
            exit 1
        else
            rotate_image.py ${image} ${tmp_image} ${rotate_degrees}
        fi
    else
        ln -s ${image} ${tmp_image}
    fi
    if tiffinfo ${tmp_image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        if od -h -j2 -N2 ${tmp_image} | head -n1 | sed 's/[0-9]*  *//' | grep -q -E '002b|2b00'
        then
            tiffcp -c none -m 0 -8 ${tmp_image} uncompressed.tif
        else
            tiffcp -c none -m 0 ${tmp_image} uncompressed.tif || tiffcp -c none -m 0 -8 ${tmp_image} uncompressed.tif
        fi
        bioformats2raw --no-hcs uncompressed.tif ${filename}.zarr
    else
        bioformats2raw --no-hcs ${tmp_image} ${filename}.zarr
    fi
    consolidate_md.py ${filename}.zarr
    """
}

process ome_zarr_metadata{
    tag "${zarr}, ${img_type}"
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

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(prefix), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem_str}-anndata.zarr"), emit: converted_anndatas, optional: true
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
    label 'build_config'
    debug verbose_log
    cache false

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(config_map), val(files), val(img_map)

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

process write_spatialdata {
    tag "${stem}"
    debug verbose_log
    
    publishDir outdir_with_version, mode: "copy"
    
    input:
    tuple val(stem), path(anndata_path), path(raw_img_path), path(label_img_path)
    
    output:
    path("${stem_str}-spatialdata.zarr")
    
    script:
    stem_str = stem.join("-")
    raw_img_str = raw_img_path ? "--raw_img_path ${raw_img_path}" : ""
    label_img_str = label_img_path ? "--label_img_path ${label_img_path}" : ""
    """
    write_spatialdata.py \
        --stem ${stem_str} \
        --anndata_path ${anndata_path} \
        ${raw_img_str} \
        ${label_img_str}
    """
}

process Generate_image {
    tag "${stem}, ${img_type}, ${file_path}"
    debug verbose_log

    publishDir outdir_with_version, mode: "copy", enabled: params.publish_generated_img

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

    warnParams()

    Process_files()

    Process_images()

    Output_to_config(
        Process_files.out.file_paths,
        Process_images.out.img_zarrs
    )
        
    if (params.write_spatialdata) {
        Output_to_spatialdata(
            Process_files.out.anndata_files,
            Process_images.out.img_tifs
        )
    }
    
}


workflow Process_files {

    warnParams()

    // Map inputs to: 
    // tuple val(stem), val(prefix), path(file), val(type), val(args)
    data_list = inputs.data.flatMap { stem, data_map ->
        data_map.data_path ?
        [
            [
                stem,
                data_map.prefix ?: "",
                file(data_map.data_path),
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
    anndata_files = route_file.out.converted_anndatas

    emit:
    files = files
    file_paths = file_paths
    anndata_files = anndata_files
}


workflow Process_images {

    warnParams()

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
            file(data_map.data_path),
            false, // keep_filename
            dataset_args[stem]?.rotate?.rotate_degrees ?: "NO_ROT" // rotate
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
            file(data_map.data_path),
            data_map.file_type,
            file(data_map.ref_img ?: "NO_REF") ,
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
                true, // keep_filename
                dataset_args[stem]?.rotate?.rotate_degrees ?: "NO_ROT" // rotate
            ]
        }
        .transpose(by: 3)
        .set {label_tifs}

    all_tifs = img_tifs.mix(label_tifs)
    all_tifs.tap{tifs}
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
    img_tifs = tifs
}


workflow Output_to_config {
    take: 
    out_file_paths
    out_img_zarrs
    
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


workflow Output_to_spatialdata {
    take: 
    anndata_files
    img_tifs
    
    main:
        img_tifs
            .map { stem, prefix, type, img, k -> 
                [stem, [type: type, img: img]]
            }
            .branch { stem, data ->
                raw: data.type == "raw"
                label: data.type == "label"
            }
        .set{tif_files}

        anndata_files
            .join(tif_files.raw, remainder: true)
            .join(tif_files.label, remainder: true)
            .map { stem, anndata, raw_tif, label_tif -> [
                stem, anndata,
                raw_tif ? raw_tif.img : [],
                label_tif ? label_tif.img : []
            ]}
            .set{data_for_sd}

        write_spatialdata(
            data_for_sd
        )
        
}

workflow {
    Full_pipeline()
}
