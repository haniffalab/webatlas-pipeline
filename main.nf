#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.title = ""
params.h5ad = ""
params.images = [
    ["image", "path-to-raw.tif"],
    ["label", "path-to-label.tif"],
]
params.max_n_worker = "30"

verbose_log = false

process h5ad_to_dict{
    tag "${h5ad}"
    echo verbose_log

    container "hamat/web-altas-data-conversion:latest"
    /*publishDir params.outdir, mode: "copy"*/
    storeDir params.outdir

    input:
        file(h5ad)

    output:
        tuple val(stem), file("*.pickle")

    script:
    stem = h5ad.baseName
    """
    h5ad_2_json.py --h5ad_file ${h5ad}
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
        --cell_sets_file cell_sets.json \
        --matrix_file matrix.json \
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

process Build_config{
    tag "config"
    echo verbose_log

    publishDir params.outdir, mode: "copy"

    output:
        file("config.json")

    script:
    """
    build_config.py --outdir ${params.outdir}
    """
}

/*
 * TODO: Build the config from from processed jsons and zarrs; Pseudo code for now
 */
process Build_config_with_md {
    tag "config"
    echo verbose_log

    publishDir params.outdir, mode: "copy"

    input:
        tuple val(stem), file(jsons)
        tuple file("*") //zarrs

    output:
        file("config.json")

    script:
    """
    build_config.py --outdir ${params.outdir} # --jsons ${jsons} --zarrs ${zarrs}
    """
}

workflow {
    h5ad_to_dict(Channel.fromPath(params.h5ad))
    dict_to_jsons(h5ad_to_dict.out)
}

workflow To_ZARR {
    channel.from(params.images)
        .map{it -> [it[0], file(it[1])]}
        .set{image_to_convert}
    image_to_zarr(image_to_convert)
}

//TODO: a one-liner to generate the json and zarr, along with the config file based on their content
workflow Full_pipeline {
    Preprocess_h5ad(Channel.fromPath(params.h5ad))
    channel.from(params.images)
        .map{it -> [it[0], file(it[1])]}
        .set{image_to_convert}
    image_to_zarr(image_to_convert)
    Build_config_with_md(Preprocess_h5ad.out, image_to_zarr.out.collect())
}

workflow Config {
    Build_config()
}

//TODO: a one-liner to generate the config file with provided jsons and zarrs
workflow Generate_config_with_processed_data {
    channel.from(params.jsons)
        .map{it -> [params.title, it]} //open to any suggestions here
        .set{jsons_with_ids}
    Build_config_with_md(jsons_with_ids, channel.fromPath(params.zarrs).collect())
}
