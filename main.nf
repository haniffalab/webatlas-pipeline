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
params.zarr_dirs = []
params.data = []
params.url = ""

params.options = []
params.layout = "minimal"
params.custom_layout = ""
params.outdir = ""

// if directly writing to s3
params.s3_keys = ["YOUR_ACCESS_KEY", "YOUR_SECRETE_KEY"]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"

params.tsv = "./template.tsv"

verbose_log = true
version = "0.0.1"

process image_to_zarr {
    tag "${image}"
    debug false

    container "openmicroscopy/bioformats2raw:0.4.0"
    storeDir params.outdir

    input:
    tuple val(stem), val(img_type), path(image)
    tuple val(accessKey), val(secretKey)
    val output_s3

    output:
    /*val out_s3, emit: s3_path*/
    path("${zarr_stem}_${img_type}.zarr"), emit: raw_zarr
    tuple val(stem), path("${zarr_stem}_${img_type}.zarr/OME/METADATA.ome.xml"), emit: ome_xml

    script:
    out_s3 = "${output_s3}/${img_type}.zarr"
    zarr_stem = file(image).baseName
    """
    #/opt/bioformats2raw/bin/bioformats2raw --output-options "s3fs_access_key=${accessKey}|s3fs_secret_key=${secretKey}|s3fs_path_style_access=true" \
        #${image} s3://${out_s3}
    /opt/bioformats2raw/bin/bioformats2raw --no-hcs ${image} ${zarr_stem}_${img_type}.zarr
    """
}

process consolidate_metadata{
    tag "${zarr}"
    /*debug verbose_log*/
    container "hamat/webatlas-zarr:${version}"

    input:
    path zarr

    script:
    """
    consolidate_md.py ${zarr}
    """
}

process any_file {
    debug verbose_log
    tag "${type}"

    container "hamat/webatlas-router"
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(type), file(file), val(args)

    output:
    file("*")

    script:
    args_strs = []
    if (args) {
        args.each { arg, value ->
            if (value instanceof Collection){
                value = value.collect { it instanceof String ? /\'/ + it.replace(" ",/\ /) + /\'/ : it }
                concat_args = value.join(',')
            }
            else
                concat_args = value
            args_strs.add("--$arg $concat_args")
        }
    }
    args_str = args_strs.join(' ')

    """
    process_${type}.py --file ${file} ${args_str}
    """
}

process route_file {
    tag "${file}"
    debug verbose_log

    container "hamat/webatlas-router"
    /*storeDir params.outdir*/
    publishDir params.outdir, mode:"copy"

    input:
    tuple val(stem), path(file), val(type), val(args)

    output:
    tuple val(stem), path("${stem}*zarr"), emit: converted_files

    script:
    args_strs = []
    if (args) {
        args.each { arg, value ->
            if (value instanceof Collection){
                value = value.collect { it instanceof String ? /\'/ + it.replace(" ",/\ /) + /\'/ : it }
                concat_args = value.join(',')
            }
            else
                concat_args = value
            args_strs.add("--$arg $concat_args")
        }
    }
    args_str = args_strs.join(' ')
    """
    router.py --type ${type} --file ${file} --stem ${stem} ${args_str}
    """
}

process Build_config{
    tag "config"
    debug verbose_log
    container "hamat/webatlas-build-config:${version}"
    publishDir params.outdir, mode: "copy"

    input:
        val(dir)
        val(title)
        val(dataset)
        val(url)
        file(zarr_dirs)
        file(files)
        val(options)
        val(layout)
        val(custom_layout)
        file(codebook)

    output:
        file("config.json")

    script:
    files = files.collect{ /\'/ + it + /\'/ }
    zarr_dirs = zarr_dirs.collect{ /\'/ + it + /\'/ }

    file_paths = files ? "--file_paths [" + files.join(',') + "]": ""
    zarr_dirs_str = zarr_dirs ? "--zarr_dirs [" + zarr_dirs.join(',') + "]" : ""
    url_str = url?.trim() ? "--url ${url}" : ""
    clayout_str = custom_layout?.trim() ? "--custom_layout \"${custom_layout}\"" : ""
    """
    build_config.py \
        --title "${title}" \
        --dataset ${dataset} \
        --files_dir ${dir} ${zarr_dirs_str} \
        --options ${options} \
        ${file_paths} ${url_str} \
        --layout ${layout} ${clayout_str} \
        --codebook ${codebook}
    """
}

process Generate_label_image {
    tag "${stem}"
    debug verbose_log

    container "generate_label:latest"
    storeDir params.outdir

    input:
        tuple val(stem), path(xml), path(h5ad), val(data_tyep), val(args)

    output:
        tuple val(stem), val("label"), file("${stem}.tif")

    script:
    """
    generate_label.py --stem "${stem}" --xml ${xml} --h5ad ${h5ad}
    """
}


Channel.fromPath(params.tsv)
    .splitCsv(header:true, sep:"\t")
    .multiMap { l ->
        images: [file(l.file).baseName, l.image_type, file(l.image_path)]
        data: [file(l.file).baseName, l.file, file(l.file).extension, l.args]
    }
    .set { data_with_md }

/*data_with_md.images.view()*/
/*data_with_md.data.view()*/

workflow {
    To_ZARR()

    Process_files()

    Generate_label_image(To_ZARR.out.all_ome_xmls.join(data_with_md.data))

    _label_to_ZARR(Generate_label_image.out)
}

workflow To_ZARR {
    if (data_with_md.images) {
        image_to_zarr(data_with_md.images, params.s3_keys, params.outdir_s3)
        consolidate_metadata(image_to_zarr.out.raw_zarr) // this will create .zmetadata in-place
        zarr_dirs = image_to_zarr.out.raw_zarr
        ome_xmls = image_to_zarr.out.ome_xml
    } else {
        zarr_dirs = []
        ome_xmls = []
    }

    emit:
        zarr_dirs = zarr_dirs
        all_ome_xmls = ome_xmls
}

workflow _label_to_ZARR {
    take: label_images

    main:
    if (label_images) {
        image_to_zarr(label_images, params.s3_keys, params.outdir_s3)
        consolidate_metadata(image_to_zarr.out.raw_zarr) // this will create .zmetadata in-place
        zarr_dirs = image_to_zarr.out.raw_zarr
    } else {
        zarr_dirs = []
    }

    emit:
        zarr_dirs = zarr_dirs
}

workflow Process_files {
    if (data_with_md.data){
        route_file(data_with_md.data)
        files = route_file.out.converted_files
    } else {
        files = []
    }

    emit:
        files = files
}

workflow Full_pipeline {

    To_ZARR()

    Process_files()

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    // Build config from files generated from Process_files
    // Ignores files in params.outdir
    Build_config(
        "''",
        params.title,
        params.dataset,
        params.url,
        To_ZARR.out.zarr_dirs.collect(),
        Process_files.out.files.collect(),
        options_str,
        params.layout,
        params.custom_layout,
        channel.fromPath(params.codebook)
    )
}

workflow Config {
    // TODO: use params.images
    if (params.zarr_dirs){
        zarr_dirs = Channel.fromPath(params.zarr_dirs).collect()
    }
    else {
        zarr_dirs = []
    }

    options_str = /"/ + new JsonBuilder(params.options).toString().replace(/"/,/\"/).replace(/'/,/\'/) + /"/

    // Build config from files in params.outdir
    Build_config(
        params.outdir,
        params.title,
        params.dataset,
        params.url,
        zarr_dirs,
        [],
        options_str,
        params.layout,
        params.custom_layout,
        channel.fromPath(params.codebook)
    )
}
