#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.title = ""
params.h5ad = ""
params.images = [
    ["image", "path-to-raw.tif"],
    ["label", "path-to-label.tif"],
]
params.factors = []
params.max_n_worker = "30"
params.dataset = ""
params.zarr_dirs = []
params.data = []

verbose_log = false

process h5ad_to_dict{
    tag "${h5ad}"
    echo verbose_log

    container "hamat/web-altas-data-conversion:latest"
    publishDir params.outdir, mode: "copy"
    // storeDir params.outdir

    input:
        file(h5ad)
        val(factors)

    output:
        tuple val(stem), file("*.pickle")

    script:
    stem = h5ad.baseName
    concat_factors = factors.join(',')
    """
    h5ad_2_json.py --h5ad_file ${h5ad} --factors ${concat_factors}
    """
}

process dict_to_jsons {
    tag "${dict}"
    echo verbose_log

    container "hamat/web-altas-data-conversion:latest"
    publishDir params.outdir, mode: "copy"

    input:
        tuple val(stem), file(dict)

    output:
        tuple val(stem), file("*.json")

    script:
    """
    dict_2_jsons.py --dict_file ${dict} \
        --cells_file cells.json \
        --cell_sets_file cell-sets.json \
        --matrix_file clusters.json \
    """
}

process image_to_zarr {
    tag "${image}"
    echo verbose_log

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
    echo true
    tag "${type}"

    conda "${type}.yaml" // conditional conda env

    input:
    tuple val(type), val(args)

    output:
    val(type)

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
    process_${type}.py ${args_str}
    """
}

process route_file {
    echo true
    tag "${type}"

    conda "global_env.yaml" // or conditional conda env?

    input:
    tuple val(type), val(args)

    output:
    val(type)

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
    router.py --type ${type} ${args_str}
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

    output:
        file("config.json")

    script:
    concat_zarr_dirs = zarr_dirs.join(',')
    """
    build_config.py --title "${title}" --dataset ${dataset} --files_dir ${dir} --zarr_dirs ${concat_zarr_dirs}
    """
}

/*
 * TODO: Build the config from from processed jsons and zarrs; Pseudo code for now
 */
process Build_config_with_md {
    tag "config"
    echo verbose_log
    containerOptions "-v ${params.outdir}:${params.outdir}"
    publishDir params.outdir, mode: "copy"

    input:
        val(dir)
        val(title)
        tuple val(stem), file(jsons)
        file(zarr_dirs)
        val(done_other_data)

    output:
        file("config.json")

    script:
    concat_zarr_dirs = zarr_dirs.join(',')
    """
    build_config.py --title "${title}" --dataset ${stem} --files_dir ${dir} --zarr_dirs ${concat_zarr_dirs}
    """
}

workflow {
    h5ad_to_dict(Channel.fromPath(params.h5ad), params.factors.collect())
    dict_to_jsons(h5ad_to_dict.out)
}

workflow To_ZARR {
    channel.from(params.images)
        .map{it -> [it[0], file(it[1])]}
        .set{image_to_convert}
    image_to_zarr(image_to_convert)
}

workflow Process_files {
    data_list = []
    params.data.each { data_type, data_map ->
        data_list.add([data_type, data_map.args])
    }
    
    route_file(Channel.from(data_list))
    emit:
        done = true
}

//TODO: a one-liner to generate the json and zarr, along with the config file based on their content
workflow Full_pipeline {
    h5ad_to_dict(Channel.fromPath(params.h5ad), params.factors.collect())
    dict_to_jsons(h5ad_to_dict.out)

    channel.from(params.images)
        .map{it -> [it[0], file(it[1])]}
        .set{image_to_convert}
    image_to_zarr(image_to_convert)
    zarr_dirs = image_to_zarr.out.collect()

    Process_files()

    Build_config_with_md(
        Channel.fromPath(params.outdir),
        params.title,
        dict_to_jsons.out,
        zarr_dirs,
        Process_files.out.done
    )
}

workflow Config {
    if (params.zarr_dirs.size > 0){
        zarr_dirs = Channel.fromPath(params.zarr_dirs).collect()
    }
    else {
        zarr_dirs = []
    }
    Build_config(
        Channel.fromPath(params.outdir),
        params.title,
        params.dataset,
        zarr_dirs
    )
}

//TODO: a one-liner to generate the config file with provided jsons and zarrs
workflow Generate_config_with_processed_data {
    channel.from(params.jsons)
        .map{it -> [params.title, it]} //open to any suggestions here
        .set{jsons_with_ids}
    Build_config_with_md(jsons_with_ids, channel.fromPath(params.zarrs).collect())
}
