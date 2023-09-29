#!/usr/bin/env/ nextflow

import groovy.json.*

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
    tuple val(dataset), path("${reindexed_label_image}")

    script:
    basename = label_image.baseName
    reindexed_label_image = "reindexed-${basename}.*"
    """
    integrate_image.py \
        --label_image ${label_image} \
        --offset ${offset} \
        --out_filename ${reindexed_label_image}
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

// @TODO
process Build_multimodal_config {
    tag "${project}"
    debug verbose_log
    cache false

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(project), val(config_map), val(datasets)

    output:
    path("${project}-multimodal-config.json")

    script:
    url_str = config_map.url?.trim() ? "--url \"${config_map.url.trim()}\"" : ""
    """
    build_config_multimodal.py \
        --project "${project}" \
        #@TODO: prepare this as json string --datasets "${datasets}" \
        --extended_features "${config_map.extend_feature_name}" \
        ${url_str} \
        --title "${config_map.title}" \
        --description "${config_map.description}"
    """
}


workflow {
    Channel.from(params.data)
        .multiMap{ it ->
            info: [it.dataset, it.obs_type ?: "cell", it.is_spatial ?: false, it.vitessce_options ?: [:]]
            raws : [it.dataset, it.raw_image] // not processed but necessary for writing config
            labels : [it.dataset, it.label_image, it.offset]
            adatas : [it.dataset, file(it.anndata), it.offset, file(it.extend_feature ?: "NO_FT")]
        }
        .set{data}
    
    // Filter null raw image
    data.raws.filter{ it[1] }.map{ [it[0], file(it[1])] }
        .set{raw_images}
    
    // Process labels filtered out nulls
    process_label(data.labels.filter { it[1] }.map{ [it[0], file(it[1]), it[2]] })

    // Reindex and concatenate anndata matrices
    // filename will be used to join output from intersection
    process_anndata(data.adatas)
    process_anndata.out
        .map{ [it[1].baseName.toString(), it[0], it[1] ] }
        .set{reid_output}

    // Perform intersection on list of anndatas
    // filename without 'intersect-' is used to match the dataset it comes from
    // as we cannot get paired dataset-file items as output
    // because output is a list of files (they're processed together)
    intersect_anndatas(process_anndata.out.collect{ it[1] })
    intersect_anndatas.out
        .flatMap()
        .map{[(it.baseName.split("intersect-")[-1]), it]}
        .set{intersect_output}
    
    // Get ome zarr metadata from raw and label images
    ome_zarr_metadata(
        raw_images.map{[
            it[0],
            "raw",
            file("${it[1].toString().replaceFirst(/\/*$/, '')}/OME/METADATA.ome.xml")
        ]}.concat(
            process_label.output.map{[
                it[0],
                "label",
                file("${it[1].toString().replaceFirst(/\/*$/, '')}/OME/METADATA.ome.xml")
            ]}
        )
    ).set{images_md}

    // Structure image data to join to corresponding dataset data
    raw_images.map{[it[0], "raw", it[1]]}
        .concat(process_label.output.map{[it[0], "label", it[1]]})
        .join(images_md, by: [0,1])
        .map { dataset, type, path, md ->
            [
                dataset, type, [path: path.name, md: new JsonSlurper().parseText(md.trim())]
            ]
        }
        .groupTuple(by: [0,1])
        .map { dataset, type, img ->
            [dataset, [type: type, img: img]]
        }
        .branch { dataset, data ->
            raw: data.type == "raw"
            label: data.type == "label"
        }
        .set{imgs}
    
    // Image data structure to be able to perform join:
    // [dataset, raw:[raw_imgs], label:[label_imgs]]
    imgs.raw
        .join(imgs.label, remainder: true)
        .map { dataset, raw_data, label_data -> [
            dataset,
            [
                raw: raw_data ? raw_data.img : [],
                label: label_data ? label_data.img : []

            ]
        ]}
        .set{img_map}
    
    // Join reindexed and intersected outputs
    // to get [dataset, intersected_file]
    reid_output
        .join(intersect_output, by: 0)
        .map{ [it[1], it[3].name] }
        .groupTuple(by: 0)
        .set{files}
    
    // @TODO
    data.info
        .join(files, by:0)
        .join(img_map, by: 0, remainder: true)
        .view()
    
    // Build_multimodal_config(
    //     data_for_config
    // )

}
