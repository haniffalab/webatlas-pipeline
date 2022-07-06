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
params.s3_keys = ["YOUR_ACCESS_KEY", "YOUR_SECRETE_KEY"]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"

params.tsv = "./template.tsv"

verbose_log = true
version = "0.0.1"

process image_to_zarr {
    tag "${image}"
    debug verbose_log

    container "openmicroscopy/bioformats2raw:0.4.0"
    storeDir params.outdir

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
    #/opt/bioformats2raw/bin/bioformats2raw --output-options "s3fs_access_key=${accessKey}|s3fs_secret_key=${secretKey}|s3fs_path_style_access=true" \
        #${image} s3://${out_s3}
    /opt/bioformats2raw/bin/bioformats2raw --no-hcs ${image} ${stem}_${img_type}.zarr
    """
}

process consolidate_metadata{
    tag "${zarr}"
    debug verbose_log
    container "hamat/webatlas-zarr:${version}"

    input:
    tuple val(stem), path(zarr)

    script:
    """
    consolidate_md.py ${zarr}
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
    /*storeDir params.outdir*/
    publishDir params.outdir, mode:"copy"

    input:
    tuple val(stem), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem}*"), emit: converted_files, optional: true

    script:
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --file_path ${file} --stem ${stem} ${args_str}
    """
}

process Build_config{
    tag "${stem}"
    debug verbose_log

    container "hamat/webatlas-build-config:${version}"
    publishDir params.outdir, mode: "copy"

    input:
        tuple val(stem), val(files), val(raster), val(label), val(raster_md), val(raw_str), val(label_md), val(label_str), val(title), val(dataset), val(url), val(options)
        val(layout)
        val(custom_layout)

    output:
        path("${stem}_config.json")

    script:
    raster = file(raster).name
    label = file(label).name
    file_paths = files.collect{ /"/ + it + /"/ }.join(",")
    url_str = url?.trim() ? "--url ${url}" : ""
    clayout_str = custom_layout?.trim() ? "--custom_layout \"${custom_layout}\"" : ""
    options_str = options ? "--options " + /"/ + new JsonBuilder(options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/ : ""

    """
    build_config.py \
        --title "${title}" \
        --dataset ${dataset} \
        --image_zarr '{"${raster}":${raster_md}, "${label}":${label_md}}' \
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
    storeDir params.outdir

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
        images: [l.title + "_" + l.dataset, l.image_type, file(l.image_path)]
        data: [
            h5ad: [l.title + "_" + l.dataset, l.h5ad, params.args.h5ad],
            molecules: [l.title + "_" + l.dataset, l.molecule, params.args.molecules]
        ]
        config_params: [l.title + "_" + l.dataset, l.title, l.dataset, l.url, l.options_str]
    }
    .set { data_with_md }

workflow ISS_pipeline {
    To_ZARR()
    /*To_ZARR.out.zarr_dirs.view()*/

    Process_files()

    To_ZARR.out.zarr_dirs
        .branch{ it ->
            raw : it[1] =~ /raw.zarr/
            label : it[1] =~ /label.zarr/
        }
        .set{zarrs}

    zarrs.raw.view()
    zarrs.label.view()

    To_ZARR.out.ome_md_json
        .branch{ it ->
            raw : it[2] ==~ /raw/
            label : it[2] ==~ /label/
        }
        .set{zarr_mds}

    /*zarr_mds.others.view()*/
    zarr_mds.raw.view()
    zarr_mds.label.view()

    Process_files.out.file_paths
            .join(zarrs.raw)
            .join(zarrs.label)
            .join(zarr_mds.raw)
            .join(zarr_mds.label)
            .join(data_with_md.config_params)
            .set{img_data_for_config}
        /*.view()*/


    Build_config(
        img_data_for_config,
        params.layout,
        params.custom_layout
    )
}

workflow To_ZARR {
    if (data_with_md.images) {
        image_to_zarr(data_with_md.images, params.s3_keys, params.outdir_s3)
        consolidate_metadata(image_to_zarr.out.raw_zarr) // this will create .zmetadata in-place
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
        consolidate_metadata(image_to_zarr.out.raw_zarr) // this will create .zmetadata in-place
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
    if (data_with_md.data){
        data_list = data_with_md.data.flatMap{ it ->
            it.collectMany{ data_type, data_map ->
                data_map[1] ? [[data_map[0], file(data_map[1]), data_type, data_map[2]]] : []
            }
        }
        route_file(data_list.unique())
        files = route_file.out.converted_files.groupTuple(by:0)
        file_paths = route_file.out.out_file_paths.map{ it -> [it[0], it[1].split('\n').flatten()].flatten() }.groupTuple(by:0)
    } else {
        files = []
        file_paths = []
    }

    emit:
        files = files
        file_paths = file_paths
}

workflow Visium_pipeline {

    To_ZARR()

    Process_files()

    h5ads = data_with_md.data.flatMap{ it ->
            it.collectMany{ data_type, data_map ->
                data_type == "h5ad" ? [[data_map[0], file(data_map[1])]] : [] }}
    Generate_label_image(To_ZARR.out.ome_md_json.join(h5ads))

    _label_to_ZARR(Generate_label_image.out)

    Process_files.out.file_paths
            .join(To_ZARR.out.zarr_dirs) //.groupTuple(by:0) if several images
            .join(_label_to_ZARR.out.label_zarr)
            .join(To_ZARR.out.ome_md_json)
            .join(_label_to_ZARR.out.ome_md_json)
            .join(data_with_md.config_params)
            .set{img_data_for_config}
        /*.view()*/


    Build_config(
        img_data_for_config,
        params.layout,
        params.custom_layout
    )
}

workflow Config {
    /*TODO: need to give an example of how to construct the args here*/
    if (params.config_files || params.zarr_md){
        Build_config(
        data_with_md.config_params,
        To_ZARR.out.zarr_dirs.collect(),
        Process_files.out.files.collect(),
        params.layout,
        params.custom_layout
        )
    }
}
