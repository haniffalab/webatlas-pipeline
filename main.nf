#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

params.outdir = "output/name"
params.h5ad = ""
params.image = ""
params.cells = ""
params.spots = ""
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
        #--genes_file ${stem}_genes.json
    """
}

process image_to_zarr {
    tag "${image}"
    echo true

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
        file(image)

    output:
        tuple val(image), file("image")

    when:
        params.image

    script:
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $image "image"
    consolidate_md.py "image/data.zarr"
    """
}

process cells_to_zarr {
    tag "${cells}"
    echo true

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
        file(cells)

    output:
        tuple val(cells), file("cells")

    when:
        params.cells

    script:
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $cells "cells"
    consolidate_md.py "cells/data.zarr"
    """
}

process spots_to_zarr {
    tag "${spots}"
    echo true

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
        file(spots)

    output:
        tuple val(spots), file("spots")

    when:
        params.spots

    script:
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $spots "spots"
    consolidate_md.py "spots/data.zarr"
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

workflow To_ZARR {
    if( params.image )
        image_to_zarr(channel.fromPath(params.image))

    if( params.cells )
        cells_to_zarr(channel.fromPath(params.cells))

    if( params.spots )
        spots_to_zarr(channel.fromPath(params.spots))                
}

workflow config {
    Build_config()
}