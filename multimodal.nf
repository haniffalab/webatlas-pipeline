#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

params.outdir = "./"
params.data_params="./templates/multimodal_template.csv"
params.data_params_delimiter=","
version="0.0.1"
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
    data = Channel.fromPath(params.data_params)
        .splitCsv(header:true, sep:params.data_params_delimiter, quote:"'")
        .multiMap{ it ->
            labels : [it.label_image, it.offset]
            adatas : [it.anndata, it.offset, file(it.features)]
        }
    
    process_label(data.labels.filter { it[0] })

    process_anndata(data.adatas)

    intersect_anndatas(process_anndata.out.collect())
}
