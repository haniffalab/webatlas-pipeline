#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

// Default params
params.max_n_worker = 30

params.title = ""
params.dataset = ""

params.tsv = "./template.tsv"
params.tsv_delimiter = "\t"

params.images = []
params.args = []
params.url = ""

params.options = []
params.layout = "minimal"
params.custom_layout = ""
params.outdir = ""
params.config_files = []

// if directly writing to s3
params.s3 = false
params.s3_keys = [
    "YOUR_ACCESS_KEY",
    "YOUR_SECRET_KEY"
]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"


verbose_log = true
version = "0.0.1"
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"

process image_to_zarr {
    tag "${image}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-image-to-zarr:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(img_type), path(image)
    tuple val(accessKey), val(secretKey)
    val output_s3

    output:
    /*val out_s3, emit: s3_path*/
    tuple val(stem), path("${stem_str}_${img_type}.zarr"), emit: img_zarr
    tuple val(stem), path("${stem_str}_${img_type}.zarr/OME/METADATA.ome.xml"), val(img_type), emit: ome_xml

    script:
    stem_str = stem.join("-")
    out_s3 = "${output_s3}/${img_type}.zarr"
    """
    #/opt/bioformats2raw/bin/bioformats2raw --output-options "s3fs_access_key=${accessKey}|s3fs_secret_key=${secretKey}|s3fs_path_style_access=true" ${image} s3://${out_s3}
    if tiffinfo ${image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        tiffcp -c none ${image} uncompressed.tif
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs uncompressed.tif ${stem_str}_${img_type}.zarr
    else
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs ${image} ${stem_str}_${img_type}.zarr
    fi
    consolidate_md.py ${stem_str}_${img_type}.zarr
    """
}

process ome_zarr_metadata{
    tag "${zarr}"
    debug verbose_log
    container "haniffalab/vitessce-pipeline-processing:${version}"

    input:
    tuple val(stem), path(zarr), val(img_type)

    output:
    tuple val(stem), stdout, val(img_type)
    /*[>tuple val(zarr), stdout<]*/

    script:
    """
    ome_zarr_metadata.py --xml_path ${zarr}
    """
}

process route_file {
    tag "${file}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem}*"), emit: converted_files, optional: true

    script:
    stem_str = stem.join("-")
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --path ${file} --stem ${stem_str} ${args_str}
    """
}

process Build_config{
    tag "${stem}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-build-config:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(files), val(raster), val(label), val(raster_md), val(raw_str), val(label_md), val(label_str), val(title), val(dataset), val(url), val(options)
    val(layout)
    val(custom_layout)

    output:
    path("${stem_str}_config.json")

    script:
    stem_str = stem.join("-")
    zarrs = [] + (raster && raster.toString() ? "\"${raster.name}\":${raster_md}" : []) + (label && label.toString() ? "\"${label.name}\":${label_md}" : [])
    zarrs_str = zarrs ? "--image_zarr '{" + zarrs.join(",") + "}'" : ""
    file_paths = files.collect{ /"/ + it + /"/ }.join(",")
    url_str = url?.trim() ? "--url ${url.replaceFirst(/\/*$/, "")}/${version}" : ""
    clayout_str = custom_layout?.trim() ? "--custom_layout \"${custom_layout}\"" : ""
    options_str = options ? "--options '" + (options instanceof String ? options : new JsonBuilder(options).toString()) + "'" : ""
    """
    build_config.py \
        --title "${title}" \
        --dataset ${dataset} \
        ${zarrs_str} \
        --file_paths '[${file_paths}]' \
        ${url_str} \
        ${options_str} \
        --layout ${layout} ${clayout_str}
    """
}

process Generate_label_image {
    tag "${stem}"
    debug verbose_log

    container "haniffalab/vitessce-pipeline-processing:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), path(file_path), path(ref_img)

    output:
    tuple val(stem), val("label"), file("${stem_str}.tif")

    script:
    stem_str = stem.join("-")
    """
    generate_label.py --stem "${stem_str}" --ome_md '${ome_md_json}' --h5ad ${h5ad}
    """
}


Channel.fromPath(params.tsv)
    .splitCsv(header:true, sep:params.tsv_delimiter, quote:"'")
    // .map { l -> tuple( "${l.title}-${l.dataset}", l ) }
    .map { l -> tuple( tuple(l.title, l.dataset), l ) }
    .branch { stem, l ->
        images: l.data_type in ["raw_image","label_image","label_image_data"]
        data: l.data_type in ["h5ad","spaceranger"]
        other: true
    }
    .set{inputs}


workflow Entry {

    // Process_files()

    Process_images()

}

workflow Process_files {
    inputs.data.view()
    if (inputs.data){
        data_list = inputs.data.flatMap { stem, data_map ->
            data_map.data_path ? [
                [
                    stem,
                    data_map.data_path,
                    data_map.data_type,
                    (data_map.args && data_map.args?.trim() ?
                        data_map.args?.trim() :
                        params.args[data_map.data_type] ?
                            params.args[data_map.data_type] :
                            []
                    )
                ]
            ] : []
        }
        route_file(data_list.unique())

        files = route_file.out.converted_files.groupTuple(by:0)
        file_paths = route_file.out.out_file_paths.map{ it -> [
                it[0],
                it[1].split('\n').flatten()
            ].flatten() }
            .groupTuple(by:0)
    } else {
        files = []
        file_paths = []
    }

    emit:
    files = files
    file_paths = file_paths
}

workflow Process_images {
    inputs.images.view()
    if (inputs.images) {
        img_tifs = inputs.images.filter { stem, data_map -> 
            data_map.data_type in ["raw_image", "label_image"] 
        }
        .map { stem, data_map -> 
            [ 
                stem,
                data_map.data_type.replace("_image",""),
                data_map.data_path
            ]
        }

        img_data = inputs.images.filter { stem, data_map ->
            data_map.data_type == "label_image_data"
        }
        .map { stem, data_map ->
            [
                stem,
                data_map.data_path,
                *(data_map.args && data_map.args?.trim() && new JsonSlurper().parseText(data_map.args).containsKey("ref_img") ? 
                    [
                        new JsonSlurper().parseText(data_map.args).ref_img,
                        new JsonSlurper().parseText(data_map.args).findAll { it.key != "ref_img" }
                    ] :
                    [
                        "NO_REF",
                        data_map.args && data_map.args?.trim() ? data_map.args : []
                    ]
                )
            ]
        }

        img_data.view()

        img_zarrs = []
        ome_md_json = []

        // Generate_label_image

        // image_to_zarr(img_tifs, params.s3_keys, params.outdir_s3)
        // img_zarrs = image_to_zarr.out.img_zarr

        // ome_zarr_metadata(image_to_zarr.out.ome_xml)
        // ome_md_json = ome_zarr_metadata.out
    } else {
        img_zarrs = []
        ome_md_json = []
    }

    emit:
    img_zarrs = img_zarrs
    ome_md_json = ome_md_json
}


workflow To_ZARR {
    if (data_with_md.images) {
        image_to_zarr(data_with_md.images, params.s3_keys, params.outdir_s3)
        zarr_dirs = image_to_zarr.out.img_zarr

        ome_zarr_metadata(image_to_zarr.out.ome_xml)
        ome_md_json = ome_zarr_metadata.out//{ [[(it[0]): new JsonSlurper().parseText(it[2].replace('\n',''))]] }
    } else {
        zarr_dirs = []
        ome_md_json = []
    }

    emit:
    zarr_dirs = zarr_dirs
    ome_md_json = ome_md_json
}

workflow _label_to_ZARR {
    take: label_images

    main:
    if (label_images) {
        image_to_zarr(label_images, params.s3_keys, params.outdir_s3)
        label_zarr = image_to_zarr.out.img_zarr
        ome_zarr_metadata(image_to_zarr.out.ome_xml)
        ome_md_json = ome_zarr_metadata.out
    } else {
        label_zarr = []
        ome_md_json = []
    }

    emit:
    label_zarr = label_zarr
    ome_md_json = ome_md_json
}


workflow scRNAseq_pipeline {
    Process_files()

    Process_files.out.file_paths
        .map{ it ->
            it + [
                null,
                null,
                null,
                null,
                null,
                null
            ] // null image paths
        }
        .join(data_with_md.config_params)
        .set{img_data_for_config}

    Build_config(
        img_data_for_config,
        params.layout,
        params.custom_layout
        )
}

workflow ISS_pipeline {
    To_ZARR()

    Process_files()

    To_ZARR.out.zarr_dirs
        .branch{ it ->
            raw : it[1] =~ /raw.zarr/
            label : it[1] =~ /label.zarr/
        }
        .set{zarrs}

    To_ZARR.out.ome_md_json
        .branch{ it ->
            raw : it[2] ==~ /raw/
            label : it[2] ==~ /label/
        }
        .set{zarr_mds}

    zarr_mds.raw.view()
    zarr_mds.label.view()

    Process_files.out.file_paths
        .join(zarrs.raw)
        .join(zarrs.label)
        .join(zarr_mds.raw)
        .join(zarr_mds.label)
        .join(data_with_md.config_params)
        .set{img_data_for_config}

    Build_config(
        img_data_for_config,
        params.layout,
        params.custom_layout
        )
}

workflow Visium_pipeline {
    To_ZARR()

    Process_files()

    h5ads = data_with_md.data.flatMap{ it ->
        it.collectMany{ data_type, data_map ->
            (data_type == "h5ad" || data_type == "spaceranger") && data_map[1] ? [[data_map[0], data_map[1]]] : []
        }
    }
    Generate_label_image(To_ZARR.out.ome_md_json.join(h5ads))

    _label_to_ZARR(Generate_label_image.out)

    Process_files.out.file_paths
        .join(To_ZARR.out.zarr_dirs) //.groupTuple(by:0) if several images
        .join(_label_to_ZARR.out.label_zarr)
        .join(To_ZARR.out.ome_md_json)
        .join(_label_to_ZARR.out.ome_md_json)
        .join(data_with_md.config_params)
        .set{img_data_for_config}

    Build_config(
        img_data_for_config,
        params.layout,
        params.custom_layout
        )
}

workflow Config {
    Build_config(
        [
            "${params.title}_${params.dataset}",
            params.files,
            new File(params.raw.zarr),
            new File(params.label.zarr),
            new JsonBuilder(params.raw.md).toString(),
            "raw",
            new JsonBuilder(params.label.md).toString(),
            "label",
            params.title,
            params.dataset,
            params.url,
            params.options
        ],
        params.layout,
        params.custom_layout
        )
}
