#!/usr/bin/env/ nextflow

// Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
import groovy.json.*

nextflow.enable.dsl=2

// Default params
params.title = ""
params.images = []
params.factors = []
params.codebook = ""
params.max_n_worker = 30
params.dataset = ""
params.zarr_md = []
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

params.tsv = "./template.tsv"

verbose_log = true
version = "0.0.1"
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"

process image_to_zarr {
    tag "${image}"
    debug verbose_log

    container "hamat/webatlas-image-to-zarr:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(img_type), path(image)
    tuple val(accessKey), val(secretKey)
    val output_s3

    output:
    /*val out_s3, emit: s3_path*/
    tuple val(stem), path("${stem}_${img_type}.zarr"), emit: raw_zarr
    tuple val(stem), path("${stem}_${img_type}.zarr/OME/METADATA.ome.xml"), val(img_type), emit: ome_xml

    script:
    out_s3 = "${output_s3}/${img_type}.zarr"
    """
    #/opt/bioformats2raw/bin/bioformats2raw --output-options "s3fs_access_key=${accessKey}|s3fs_secret_key=${secretKey}|s3fs_path_style_access=true" ${image} s3://${out_s3}
    if tiffinfo ${image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        tiffcp -c none ${image} uncompressed.tif
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs uncompressed.tif ${stem}_${img_type}.zarr
    else
        /opt/bioformats2raw/bin/bioformats2raw --no-hcs ${image} ${stem}_${img_type}.zarr
    fi
    consolidate_md.py ${stem}_${img_type}.zarr
    """
}

process ome_zarr_metadata{
    tag "${zarr}"
    debug verbose_log
    container "hamat/webatlas-ome-zarr-metadata:${version}"

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

    container "hamat/webatlas-router:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem}*"), emit: converted_files, optional: true

    script:
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --path ${file} --stem ${stem} ${args_str}
    """
}

process Build_config{
    tag "${stem}"
    debug verbose_log

    container "hamat/webatlas-build-config:${version}"
    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(files), val(raster), val(label), val(raster_md), val(raw_str), val(label_md), val(label_str), val(title), val(dataset), val(url), val(options)
    val(layout)
    val(custom_layout)

    output:
    path("${stem}_config.json")

    script:
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

    container "hamat/webatlas-generate-label:${version}"
    publishDir outdir_with_version, mode:"copy"

    input:
    tuple val(stem), val(ome_md_json), val(img_type), path(h5ad)

    output:
    tuple val(stem), val("label"), file("${stem}.tif")

    script:
    """
    generate_label.py --stem "${stem}" --ome_md '${ome_md_json}' --h5ad ${h5ad}
    """
}


Channel.fromPath(params.tsv)
    .splitCsv(header:true, sep:"\t")
    .multiMap { l ->
        images: [
            "${l.title}_${l.dataset}",
            l.image_type,
            l.image_path
        ]
        data: [
            h5ad: [
                "${l.title}_${l.dataset}",
                l.h5ad,
                l.args && l.args.h5ad?.trim() ? l.args.h5ad?.trim() : params.args.h5ad
            ],
            molecules: [
                "${l.title}_${l.dataset}",
                l.molecules,
                l.args && l.args.molecules?.trim() ? l.args.molecules?.trim() : params.args.molecules
            ],
            spaceranger: [
                "${l.title}_${l.dataset}",
                l.spaceranger,
                (l.args && l.args.spaceranger?.trim() ? l.args.spaceranger?.trim() : params.args.spaceranger ? params.args.spaceranger : []) +
                (l.args && l.args.h5ad?.trim() ? l.args.h5ad?.trim() : params.args.h5ad ? params.args.h5ad : [])
            ]
        ]
        config_params: [
            "${l.title}_${l.dataset}",
            l.title,
            l.dataset,
            l.url,
            l.options?.trim() ? l.options : params.options
        ]
    }
    .set { data_with_md }


workflow To_ZARR {
    if (data_with_md.images) {
        image_to_zarr(data_with_md.images, params.s3_keys, params.outdir_s3)
        zarr_dirs = image_to_zarr.out.raw_zarr

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
        label_zarr = image_to_zarr.out.raw_zarr
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

workflow Process_files {
    data_with_md.data.view()
    if (data_with_md.data){
        data_list = data_with_md.data.flatMap{ it ->
            it.collectMany{ data_type, data_map ->
                data_map[1] ? [
                    [
                        data_map[0],
                        data_map[1],
                        data_type,
                        data_map[2]
                    ]
                ] : []
            }
        }
        route_file(data_list.unique())
        files = route_file.out.converted_files.groupTuple(by:0)
        file_paths = route_file.out.out_file_paths.map{ it -> [
                it[0],
                it[1].split('\n').flatten()
            ].flatten() }.groupTuple(by:0)
    } else {
        files = []
        file_paths = []
    }

    emit:
    files = files
    file_paths = file_paths
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
