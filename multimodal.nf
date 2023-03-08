#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

params.outdir = "./"
params.data_params="./templates/multimodal_template.csv"
params.data_params_delimiter=","
version="0.0.1"
verbose_log=true
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"


process process_label {
    tag "${type}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"


    input:
    tuple val(type), val(offset), path(label_tif)


    output:
    tuple val(type), path("${reindexed_tif}")

    script:
    stem = label_tif.baseName
    reindexed_tif = "${stem}_${type}_reindexed.tif"
    """
    integrate_image.py --label_image $label_tif --offset ${offset} --out_filename ${reindexed_tif}
    """
}


process process_h5ads {
    tag "${type}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"


    input:
    tuple val(type), val(offset), path(label_tif)


    output:
    tuple val(type), path("${reindexed_tif}")

    script:
    """
    """
}


process intersection {
    tag "${type}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"


    input:
    tuple val(type), val(offset), path(label_tif)


    output:
    tuple val(type), path("${reindexed_tif}")

    script:
    """
    """
}


workflow {
    datas = Channel.fromPath(params.data_params)
        .splitCsv(header:true, sep:params.data_params_delimiter, quote:"'")
        .multiMap{ it ->
            labels : [it.type, it.offset, it.label_tif]
            h5ads : [it.type, it.offset, it.h5ad]
        }
    /*datas.labels.view()*/
    /*datas.h5ads.view()*/
    process_label(datas.labels)
    /*process_h5ads(datas.h5ads)*/
    /*intersection(process_h5ads.out.collect())*/
}
