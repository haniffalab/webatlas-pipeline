#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

// Default params
params.title = ""
params.images = [
    ["image", "path-to-raw.tif"],
    ["label", "path-to-label.tif"],
]
params.factors = []
params.max_n_worker = "30"
params.dataset = ""
params.zarr_dirs = []
params.data = []
params.options = []


verbose_log = true


process image_to_zarr {
    tag "${image}"
    echo false

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(img_type), file(image)

    output:
    file(img_type)

    script:
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $image "${img_type}"
    consolidate_md.py "${img_type}/data.zarr"
    """
}

process any_file {
    echo verbose_log
    tag "${type}"

    conda "global_env.yaml" // or conditional conda env
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(type), file(file), val(args)

    output:
    file("*")

    script:
    args_strs = []
    if (args) {
        args.each { arg, value ->
            if (value instanceof Collection){
                value = value.collect { it instanceof String ? /\'/ + it.replace(" ",/\ /) + /\'/ : it }
                concat_args = value.join(',')
            }
            else
                concat_args = value
            args_strs.add("--$arg $concat_args")
        }
    }
    args_str = args_strs.join(' ')

    """
    process_${type}.py --file ${file} ${args_str}
    """
}

process route_file {
    echo verbose_log
    tag "${type}"

    conda "global_env.yaml" // or conditional conda env
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(type), file(file), val(args)

    output:
    tuple val(type), file("*")

    script:
    args_strs = []
    if (args) {
        args.each { arg, value ->
            if (value instanceof Collection){
                value = value.collect { it instanceof String ? /\'/ + it.replace(" ",/\ /) + /\'/ : it }
                concat_args = value.join(',')
            }
            else
                concat_args = value
            args_strs.add("--$arg $concat_args")
        }
    }
    args_str = args_strs.join(' ')

    """
    router.py --type ${type} --file ${file} ${args_str}
    """
}

process Build_config{
    tag "config"
    echo verbose_log
    containerOptions "-v ${params.outdir}:${params.outdir}"
    publishDir params.outdir, mode: "copy"

    input:
        val(dir)
        val(title)
        val(dataset)
        file(zarr_dirs)
        val(files)
        val(options)

    output:
        file("config.json")

    script:
    concat_zarr_dirs = zarr_dirs.join(',')
    """
    build_config.py --title "${title}" --dataset ${dataset} --files_dir ${dir} --zarr_dirs ${concat_zarr_dirs} --options ${options}
    """
}

workflow {
    Process_files()
    // Process_files.out.files.toList().view()
}

workflow To_ZARR {
    channel.from(params.images)
        .map{it -> [it[0], file(it[1])]}
        .set{image_to_convert}
    image_to_zarr(image_to_convert)

    emit:
        zarr_dirs = image_to_zarr.out.collect()
}

workflow Process_files {
    data_list = []
    params.data.each { data_type, data_map ->
        data_list.add([data_type, file(data_map.file), data_map.args])
    }
    
    route_file(Channel.from(data_list))

    emit:
        files = route_file.out
}

workflow Full_pipeline {
    To_ZARR()

    Process_files()

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    Build_config(
        Channel.fromPath(params.outdir),
        params.title,
        params.dataset,
        To_ZARR.out.zarr_dirs,
        Process_files.out.files.toList(),
        options_str
    )
}

workflow Config {
    if (params.zarr_dirs.size > 0){
        zarr_dirs = Channel.fromPath(params.zarr_dirs).collect()
    }
    else {
        zarr_dirs = []
    }

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    Build_config(
        Channel.fromPath(params.outdir),
        params.title,
        params.dataset,
        zarr_dirs,
        [],
        options_str
    )
}
