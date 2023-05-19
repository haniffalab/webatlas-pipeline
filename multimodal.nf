#!/usr/bin/env/ nextflow

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
    tuple path(label_image), val(offset)

    output:
    path("${reindexed_tif}")

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
    tuple path(anndata), val(offset), path(features)

    output:
    path("*")

    script:
    """
    integrate_anndata.py reindex_and_concat \
        --path ${anndata} \
        --offset ${offset} \
        --features ${features}
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
            labels : [it.label_image, it.offset]
            adatas : [it.anndata, it.offset, file(it.features)]
        }
        .set{data}
    
    process_label(data.labels.filter { it[0] })

    process_anndata(data.adatas)

    intersect_anndatas(process_anndata.out.collect())
}
