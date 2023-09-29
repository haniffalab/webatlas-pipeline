#!/usr/bin/env/ nextflow

include { ome_zarr_metadata } from "./main.nf"

nextflow.enable.dsl=2

params.outdir = ""

version="0.3.2"
verbose_log=true
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"


process process_label {
    tag "${label_image}"
    debug verbose_log

    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(dataset), path(label_image), val(offset)

    output:
    tuple val(dataset), path("${reindexed_tif}")

    script:
    basename = label_image.baseName
    reindexed_tif = "reindexed-${basename}.tif"
    """
    integrate_image.py \
        --label_image ${label_image} \
        --offset ${offset} \
        --out_filename ${reindexed_tif}
    """
}


process process_anndata {
    tag "${anndata}"
    debug verbose_log

    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(dataset), path(anndata), val(offset), path(features)

    output:
    tuple val(dataset), path("*")

    script:
    features_str = features.name != "NO_FT" ? "--features ${features}" : ""
    """
    integrate_anndata.py reindex_and_concat \
        --path ${anndata} \
        --offset ${offset} \
        ${features_str}
    """
}


process intersect_anndatas {
    debug verbose_log

    publishDir outdir_with_version, mode:"copy"

    input:
    path(anndatas)

    output:
    path("*")

    script:
    """
    integrate_anndata.py intersect_features --paths ${anndatas}
    """
}


workflow {
    Channel.from(params.data)
        .multiMap{ it ->
            raws : [it.dataset, it.raw_image] // not processed but necessary for writing config
            labels : [it.dataset, it.label_image, it.offset]
            adatas : [it.dataset, file(it.anndata), it.offset, file(it.extend_feature ?: "NO_FT")]
        }
        .set{data}
    
    data.raws.filter{ it[1] }.map{ [it[0], file(it[1])] }
        .set{raw_images}
        
    process_label(data.labels.filter { it[1] }.map{ [it[0], file(it[1]), it[2]] })

    process_anndata(data.adatas)
    process_anndata.out
        .map{ [it[1].baseName.toString(), it[0], it[1] ] }
        .set{reid_output}

    intersect_anndatas(process_anndata.out.collect{ it[1] })
    intersect_anndatas.out
        .flatMap()
        .map{[(it.baseName.split("intersect-")[-1]), it]}
        .set{intersect_output}
    
    reid_output.join(intersect_output, by: 0)
        .map{ [it[1], it[3]] }
        .join(raw_images, by: 0, remainder: true)
        .join(process_label.output, by: 0, remainder: true)
        .view()

    // @TODO: ome_zarr_medatada()
}
