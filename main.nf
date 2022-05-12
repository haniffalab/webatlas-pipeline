#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

// Default params
params.title = ""
params.images = []
params.factors = []
params.max_n_worker = "30"
params.dataset = ""
params.zarr_dirs = []
params.data = []
params.url = ""
params.options = []


verbose_log = true


process image_to_zarr {
    tag "${image}"
    debug false

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
    debug verbose_log
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
    debug verbose_log
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
    router.py --type ${type} --file ${file} ${args_str}
    """
}

process Build_config{
    tag "config"
    debug verbose_log
    containerOptions "-v ${params.outdir}:${params.outdir}"
    publishDir params.outdir, mode: "copy"

    input:
        val(dir)
        val(title)
        val(dataset)
        val(url)
        file(zarr_dirs)
        file(files)
        val(options)

    output:
        file("config.json")

    script:
    files = files.collect{ /\'/ + it + /\'/ }
    zarr_dirs = zarr_dirs.collect{ /\'/ + it + /\'/ }

    file_paths = files ? "--file_paths [" + files.join(',') + "]": ""
    zarr_dirs_str = zarr_dirs ? "--zarr_dirs [" + zarr_dirs.join(',') + "]" : ""
    url_str = url?.trim() ? "--url ${url}" : ""
    """
    build_config.py --title "${title}" --dataset ${dataset} --files_dir ${dir} ${zarr_dirs_str} --options ${options} ${file_paths} ${url_str}
    """
}

workflow {
    Process_files()
    // Process_files.out.files.toList().view()
}

workflow To_ZARR {
    if (params.images) {
        channel.from(params.images)
            .map{it -> [it[0], file(it[1])]}
            .set{image_to_convert}
        image_to_zarr(image_to_convert)
        zarr_dirs = image_to_zarr.out.collect()
    }
    else
        zarr_dirs = []
    
    emit:
        zarr_dirs = zarr_dirs
}

workflow Process_files {
    if (params.data){
        data_list = []
        params.data.each { data_type, data_map ->
            data_list.add([data_type, file(data_map.file), data_map.args])
        }
        route_file(Channel.from(data_list))
        files = route_file.out.collect()
    }
    else
        files = []
    
    emit:
        files = files
}

workflow Full_pipeline {
    To_ZARR()

    Process_files()

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    // Build config from files generated from Process_files
    // Ignores files in params.outdir
    Build_config(
        "''",
        params.title,
        params.dataset,
        params.url,
        To_ZARR.out.zarr_dirs,
        Process_files.out.files,
        options_str
    )
}

workflow Config {
    // TODO: use params.images
    if (params.zarr_dirs){
        zarr_dirs = Channel.fromPath(params.zarr_dirs).collect()
    }
    else {
        zarr_dirs = []
    }

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    // Build config from files in params.outdir
    Build_config(
        params.outdir,
        params.title,
        params.dataset,
        params.url,
        zarr_dirs,
        [],
        options_str
    )
}
