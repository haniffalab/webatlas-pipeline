#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.title = ""
params.h5ad = ""
params.image = ""
params.labels = ""
params.max_n_worker = "30"

process Preprocess_h5ad{
    tag "${h5ad}"
    echo true

    container "hamat/web-altas-data-conversion:latest"
    publishDir params.outdir, mode: "copy"

    input:
        file(h5ad)

    output:
        tuple val(stem), file("*.json")

    script:
    stem = h5ad.baseName
    """
    h5ad_2_json.py --h5ad_file ${h5ad} \
        --cells_file cells.json \
        --cell_sets_file cell_sets.json \
        --matrix_file matrix.json \
    """
}

process To_zarr {
    tag "${params.zarr_type}"
    echo true

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
        file(source)

    output:
        tuple val(source), file("${params.zarr_type}")

    script:
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $source "${params.zarr_type}"
    consolidate_md.py "${params.zarr_type}/data.zarr"
    """
}

process Build_config{
    tag "config"
    echo true

    publishDir params.outdir, mode: "copy"

    output:
        file("config.json")

    script:
    """
    build_config.py --outdir ${params.outdir}
    """
}

workflow {
    Preprocess_h5ad(Channel.fromPath(params.h5ad))
}

workflow Image {
    if( params.image )
        params.zarr_type = "image"
        To_zarr(channel.fromPath(params.image))   
             
}

workflow Labels {
    if( params.labels )
        params.zarr_type = "labels"
        To_zarr(channel.fromPath(params.labels))   
             
}

workflow Config {
    Build_config()
}