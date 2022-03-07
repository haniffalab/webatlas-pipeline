#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>

nextflow.enable.dsl=2

/*params.h5ad = "/nfs/team283_imaging/playground_Tong/Webaltas_data/h5ad_with_cell_contours/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet.h5ad"*/
params.h5ad = "./data-files/visium.h5ad"
params.outdir = "./test"
params.img = "/home/ubuntu/Documents/planer-nf/test/out_opt_flow_registered_DAPI_flow_lab.tif"

params.max_n_worker = 30


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
        --cells_file ${stem}_cells.json \
        --cell_sets_file ${stem}_cell_sets.json \
        --matrix_file ${stem}_matrix.json \
        #--genes_file ${stem}_genes.json
    """
}

process to_zarr {
    tag "${img}"
    echo true

    conda "zarr_convert.yaml"
    publishDir params.outdir, mode: "copy"

    input:
    file(img)

    output:
    tuple val(stem), file("$stem")

    script:
    stem = img.baseName
    """
    bioformats2raw --max_workers ${params.max_n_worker} --resolutions 7 --file_type zarr $img "${stem}"
    consolidate_md.py "${stem}/data.zarr"
    """
}

workflow {
    Preprocess_h5ad(Channel.fromPath(params.h5ad))
}

workflow To_ZARR {
    to_zarr(Channel.fromPath(params.img))
}

