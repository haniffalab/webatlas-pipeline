#!/usr/bin/env/ nextflow

import groovy.json.*

nextflow.enable.dsl=2

// This should've been false by default
verbose_log = true
version = "0.5.3"

//////////////////////////////////////////////////////

// Default params
params.max_n_worker = 30

params.outdir = ""
params.args = [:]
params.projects = []
params.write_spatialdata = false
params.publish_generated_img = false

params.vitessce_options = [:]
params.layout = "minimal"
params.custom_layout = ""

params.vitessce_config_map = [
    url: "http://localhost:3000/",
    options: params.vitessce_options,
    layout: params.layout,
    custom_layout: params.custom_layout,
    title: "",
    description: ""
]

// if directly writing to s3
params.s3 = false
params.s3_keys = [
    "YOUR_ACCESS_KEY",
    "YOUR_SECRET_KEY"
]
params.outdir_s3 = "cog.sanger.ac.uk/webatlas/"

//////////////////////////////////////////////////////

// Valid types
data_types = ["h5ad","spaceranger","xenium","merscope","molecules"]
image_types = ["raw_image","label_image","raw_image_data","label_image_data"]

// Valid vitessce params
vitessce_params = ["title","description","url","vitessce_options","layout","custom_layout"]

// Add version to outdir as subdir, handling trailing slashes
outdir_with_version = "${params.outdir.replaceFirst(/\/*$/, "")}\/${version}"

//////////////////////////////////////////////////////

// Organize input params into channels to send to processes
// Workflows will further structure channels depending on data type
// This block of code outside a workflow is executed when running any workflow of the pipeline, even when running other .nf files
// so it needs params.projects = [] defined above, otherwise running another .nf file would hang while waiting for params.projects (?)

// Expected input from params file
// projects:
//   - project: "project_name"
//     args:
//       h5ad:
//         var_index: "SYMBOL"
//     datasets:
//       - dataset: "dataset_name"
//         args:
//           h5ad:
//             var_index: "SYMBOL"
//         data:
//           data_type: "h5ad"
//           data_path: "path/to/data"

// Expected output
// datasets.data channel with values per dataset as [tuple(project, dataset), data]
// datasets.config_map channel with values per dataset as [tuple(project, dataset), config_map] (for vitessce config)

// Get params inputs to per dataset channel structure
Channel.from(params.projects) // Create a channel from the projects
    .map { p -> [p.project, p.datasets] } // Map each project to tuples [project, project_datasets(array)]
    .transpose() // Flatten the tuples to [project, project_dataset1], [project, project_dataset2], ...
    .multiMap { // Create two separate channels: data and config_map. This to have vitessce params available until vitessce config generation and not carry over with each data processing
        project, dataset -> 
            data: [ tuple(project, dataset.dataset), dataset.data ] // tuple project and dataset to be used to join afterwards
            config_map: [
                tuple(project, dataset.dataset), // tuple project and dataset to be used to join afterwards
                params.vitessce_config_map + dataset.subMap(vitessce_params)  // merge vitessce_config_map in root of params file with vitessce_params defined within dataset
            ]
    }
    .set {datasets}

//Handling the raw_img_path defined in args for the data item. If data type is spaceranger - if raw_img_path exists store this, if not infer using data_path (assumes only 1 tif in data_path).
datasets.data
    .transpose(by:1)
    .filter { it[1].data_type in ['spaceranger', 'xenium'] }  // Ensure data_type is 'spaceranger'
    .map { item -> 
        def metadata = item[0]
        def dataInfo = item[1]
        def rawImgPath = dataInfo.args?.raw_img_path ?: dataInfo.data_path  // Use raw_img_path if present, otherwise use data_path

        return [ metadata, [ data_type: 'raw_image', data_path: rawImgPath ] ]
    }
    .set{raw_images_for_spaceranger}


// Expected output
// inputs.data channel with values per data as [tuple(project, dataset), *data] where data_type is in data_types
// inputs.images channel with values per data as [tuple(project, dataset), *data] where data_type is in image_types
// inputs.other channel with values per data as [tuple(project, dataset), *data] where data_type is not in data_types or image_types

// Get dataset channel to per data channel structure
datasets.data
    .transpose(by:1) // Flatten the datasets data to [tuple(project, dataset), data1], [tuple(project, dataset), data2], ...
    .branch{ stem, d -> // Split data into three channels depending on data_type: data, images and other
        data: d.data_type in data_types
        images: d.data_type in image_types
        other: true
    }
    .set{inputs}

//Ensure images used are both the ones defined in the images branch and the ones identified as data args.
all_images = inputs.images.mix(raw_images_for_spaceranger)

// If anything is in inputs.other, display error message about unrecognized data_type
inputs.other
    .collect { stem, d -> d.data_type } // Collect all data_types from inputs.other to show one error message per unique data_type
    .view{ "Unrecognized data_type(s) ${it.unique()}" }

//////////////////

// Intermediate data types
// These data types are firstly converted to h5ad before further processing
// thus this map is used to ensure the h5ad args apply to them
interm_dt = [
    spaceranger: ["h5ad"],
    xenium: ["h5ad"],
    merscope: ["h5ad"]
]

// Initialize project and dataset args as empty maps
project_args = [:]
dataset_args = [:]

// Fill project and dataset args with the provided args in the params file
// If no args are provided, default to empty map
// This allows for easy access to args in the workflow as project_args[project]
// and dataset_args[[project, dataset]]
params.projects.each{ p ->
    project_args[p.project] = p.args ?: [:]
    p.datasets.each{ d ->
        dataset_args[[p.project, d.dataset]] = d.args ?: [:]
    }
}

// Helper function to access nested map values. just from how groovy maps work
def getSubMapValues (m, keys) {
    m.subMap(keys).values().sum() ?: [:]
}

// Function used within workflows
// Merge args from root params file args, then project, then dataset, then data
// Thus data args take higher priority
// Use interm_dt to ensure h5ad args also apply to data_types that are first converted to h5ad as above
def mergeArgs (stem, data_type, args) {
    getSubMapValues(params.args, [data_type, *interm_dt[data_type]]) + 
    getSubMapValues(project_args[stem[0]], [data_type, *interm_dt[data_type]]) + 
    getSubMapValues(dataset_args[stem], [data_type, *interm_dt[data_type]]) + 
    (args ?: [:])
}

//////////////////////////////////////////////////////

// Give warning if -params-file flag is used with no file provided
// Checks for flag because params can be provided directly on the command line
def warnParams () {
    if (!workflow.commandLine.contains("-params-file")){
        log.warn "No -params-file provided"
    }
}

//////////////////////////////////////////////////////

// Use bioformats2raw to convert images to zarr
process image_to_zarr {
    tag "${image}"
    debug verbose_log

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(prefix), val(img_type), path(image), val(keep_filename)

    output:
    tuple val(stem), val(img_type), path("${filename}.zarr"), emit: img_zarr
    tuple val(stem), val(img_type), path("${filename}.zarr/OME/METADATA.ome.xml"), emit: ome_xml

    script:
    filename = keep_filename ? image.baseName : ([*stem, prefix, img_type] - null - "").join("-") // `[*stem, prefix, img_type] - null - ""` will remove nulls and empty strings from array before joining, e.g. for when no prefix is provided
    """
    if tiffinfo ${image} | grep "Compression Scheme:" | grep -wq "JPEG"
    then
        if od -h -j2 -N2 ${image} | head -n1 | sed 's/[0-9]*  *//' | grep -q -E '002b|2b00'
        then
            tiffcp -c none -m 0 -8 ${image} uncompressed.tif
        else
            tiffcp -c none -m 0 ${image} uncompressed.tif || tiffcp -c none -m 0 -8 ${image} uncompressed.tif
        fi
        bioformats2raw --no-hcs uncompressed.tif ${filename}.zarr
    else
        bioformats2raw --no-hcs ${image} ${filename}.zarr
    fi
    consolidate_md.py ${filename}.zarr
    """
}

// Read ome metadata from zarr and output to stdout
process ome_zarr_metadata{
    tag "${zarr}, ${img_type}"
    debug verbose_log

    input:
    tuple val(stem), val(img_type), path(zarr)

    output:
    tuple val(stem), val(img_type), stdout

    script:
    """
    ome_zarr_metadata.py --xml_path ${zarr}
    """
}

// Call router.py to convert data files to anndata-zarr
// args are passed as a json string for fire to read correctly
// Note args and kwargs are usually passed from router to processing functions to make it easily extensible
// thus args are named distinctively in the different functions
process route_file {
    tag "${type}, ${file}"
    debug verbose_log
    cache "lenient"

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(prefix), path(file), val(type), val(args)

    output:
    tuple val(stem), stdout, emit: out_file_paths
    tuple val(stem), path("${stem_str}-anndata.zarr"), emit: converted_anndatas, optional: true // only molecules don't output anndata-zarr, they output json, thus the anndata-zarr is optional. But probably now supported in vitessce to get molecules from anndata-zarr(?)
    tuple val(stem), path("${stem_str}*"), emit: converted_files, optional: true // both anndata-zarrs and other files (molecules json)
    tuple val(stem), path("tmp-${stem_str}*"), emit: extra_files, optional: true // for intermediary h5ad files from non-h5ad inputs

    script:
    stem_str = ([*stem, prefix] - null - "").join("-")
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}"
    """
    router.py --file_type ${type} --path ${file} --stem ${stem_str} ${args_str}
    """
}

// Call build_config.py to generate a vitessce config file
// only the string values of file (non-image) paths and image paths are provided and not the actual files as they are not needed
process Build_config {
    tag "${stem}"
    label 'build_config'
    debug verbose_log
    cache false

    publishDir outdir_with_version, mode: "copy"

    input:
    tuple val(stem), val(config_map), val(files), val(img_map)

    output:
    path("${stem_str}-config.json")

    script:
    stem_str = stem.join("-")
    file_paths = files.collect{ /"/ + it + /"/ }.join(",")
    imgs_str = img_map ? "--images '" + new JsonBuilder(img_map).toString() + "'" : ""
    url_str = config_map.url?.trim() ? "--url \"${config_map.url.trim()}\"" : ""
    options_str = config_map.options ? "--options '" + (config_map.options instanceof String ? options : new JsonBuilder(config_map.options).toString()) + "'" : ""
    clayout_str = config_map.custom_layout?.trim() ? "--custom_layout \"${config_map.custom_layout}\"" : ""
    """
    build_config.py \
        --project "${stem[0]}" \
        --dataset "${stem[1]}" \
        --file_paths '[${file_paths}]' \
        ${imgs_str} \
        ${url_str} \
        ${options_str} \
        --layout "${config_map.layout}" ${clayout_str} \
        --title "${config_map.title}" \
        --description "${config_map.description}"
    """
}

// Optional spatialdata output
// Though format was still in development at the time of writing so it is probably outdated
process write_spatialdata {
    tag "${stem}"
    debug verbose_log
    
    publishDir outdir_with_version, mode: "copy"
    
    input:
    tuple val(stem), path(anndata_path), path(raw_img_path), path(label_img_path)
    
    output:
    path("${stem_str}-spatialdata.zarr")
    
    script:
    stem_str = stem.join("-")
    raw_img_str = raw_img_path ? "--raw_img_path ${raw_img_path}" : ""
    label_img_str = label_img_path ? "--label_img_path ${label_img_path}" : ""
    """
    write_spatialdata.py \
        --stem ${stem_str} \
        --anndata_path ${anndata_path} \
        ${raw_img_str} \
        ${label_img_str}
    """
}

// Mostly for label images, generated from coordinates
// but also for merscope raw image as datasets contain individual tif images per channel so this process concatenates them
// For label images, requires either a reference image or `shape` tuple in args
// as some coordinate data for labels don't include any image shape data (e.g. xenium)
// arguably, the label image could just be extended to the maximum coordinates without needing the exact shape of the raw image
// that could be the default when no ref_img or shape is provided
process Generate_image {
    tag "${stem}, ${img_type}, ${file_path}"
    debug verbose_log

    publishDir outdir_with_version, mode: "copy", enabled: params.publish_generated_img

    input:
    tuple val(stem), val(prefix), val(img_type), path(file_path), val(file_type), path(ref_img), val(args)

    output:
    tuple val(stem), val(prefix), val(img_type), path("${stem_str}*.tif")

    script:
    stem_str = ([*stem, prefix] - null - "").join("-")
    ref_img_str = ref_img.name != "NO_REF" ? "--ref_img ${ref_img}" : ""
    args_str = args ? "--args '" + new JsonBuilder(args).toString() + "'" : "--args {}" // optional input file pattern
    """
    generate_image.py \
        --stem ${stem_str} \
        --img_type ${img_type} \
        --file_type ${file_type} \
        --file_path ${file_path} \
        ${ref_img_str} ${args_str}
    """
}

//////////////////////////////////////////////////////

// Full pipeline runs both Process_files and Process_images
// They can run in parallel as they don't depend on each other
// however, that means some processes are duplicated:
// when generating label images from spaceranger, obs can be subsetted and/or spot size can be multiplied (used in skin project where spots were merged)
// therefore to get the matching coordinates of labels the spaceranger data was first read into anndata before filtering/scaling
// otherwise, reading only the coordinates data to generate the label image without the obs metadata, we wouldn't know which coordinates match which filtered/scaled spots
// thus, this runs within Process_images the same spaceranger_to_anndata function that is run in Process_files
workflow Full_pipeline {

    warnParams()

    Process_files()

    Process_images()

    // Generating vitessce config file needs to wait for output files and images paths and metadata
    Output_to_config(
        Process_files.out.file_paths,
        Process_images.out.img_zarrs
    )
        
    if (params.write_spatialdata) {
        Output_to_spatialdata(
            Process_files.out.anndata_files,
            Process_images.out.img_tifs
        )
    }
    
}


workflow Process_files {

    warnParams()

    // Map inputs to: 
    // tuple val(stem), val(prefix), path(file), val(type), val(args)
    data_list = inputs.data.flatMap { stem, data_map ->
        data_map.data_path ?
        [
            [
                stem,
                data_map.prefix ?: "",
                file(data_map.data_path),
                data_map.data_type,
                mergeArgs(stem, data_map.data_type, data_map.args) // merge args from params, project, dataset and data
            ]
        ] : [:]
    }

    // call route_file on each data_list item
    route_file(data_list)

    files = route_file.out.converted_files // anndata-zarr and other (non tmp) files
        .map { stem, paths ->
            [ stem, [paths].flatten() ]
        } // make [ stem, [path1, path2, ...] ] <- item per route_file call, so anndata-zarr would be different item from a molecules json
        .transpose(by: 1) // make an item for each path like [ stem, path ]
        .groupTuple(by:0) // group by stem like [ stem, [path1, path2] ] across all route_file.out.converted_files, so an anndata-zarr and molecules json would be in the same item
    file_paths = files.map { stem, it -> // same as files but with paths as strings for vitessce config
        [ stem, it.name ]
    }
    anndata_files = route_file.out.converted_anndatas

    emit:
    files = files
    file_paths = file_paths
    anndata_files = anndata_files
}


workflow Process_images {

    warnParams()

    // Map tif inputs to:
    // tuple val(stem), val(prefix), val(img_type), path(image)
    img_tifs = all_images.filter { stem, data_map ->
        data_map.data_type in ["raw_image", "label_image"] // only already existing images
    }

    //Perform checks for the raw/label_image
    .map { stem, data_map -> // structure as image_to_zarr input
        def filePath = file(data_map.data_path)

        if (!filePath.exists()) {
            log.error "File not found: ${filePath}"
            throw new IllegalStateException("File not found: ${filePath}")
        }

        if (filePath.isDirectory()) {
            // Get all .tiff/.tif files in the directory
            def tiffFiles = filePath.listFiles().findAll { it.name.toLowerCase().endsWith('.tiff') || it.name.toLowerCase().endsWith('.tif') }

            if (tiffFiles.isEmpty()) {
                log.error "No .tiff files found in directory: ${filePath}"
                throw new IllegalStateException("No .tiff files found in directory: ${filePath}")
            }
            if (tiffFiles.size() > 1) {
                log.error "Multiple .tiff files found in directory: ${filePath}. Expected only one."
                throw new IllegalStateException("Multiple .tiff files found in directory: ${filePath}. Expected only one.")
            }

            // Use the single .tiff file found
            filePath = tiffFiles.first()
        }

        // Ensure the selected path is a .tiff file
        if (!filePath.name.toLowerCase().endsWith('.tiff') && !filePath.name.toLowerCase().endsWith('.tif')) {
            log.error "Invalid file format: ${filePath}. Expected .tiff or .tif"
            throw new IllegalStateException("Invalid file format: ${filePath}. Expected .tiff or .tif")
        }

        // If checks pass, return the structured output
        return [
            stem,
            data_map.prefix,
            data_map.data_type.replace("_image", ""),  // make type only `raw` or `label`
            filePath,
            false // keep_filename, output name will be generated in process as stem-prefix-type
        ]
    }

    // Map raw/label data inputs to:
    // tuple val(stem), val(prefix), val(img_type), path(file_path), val(file_type), path(ref_img), val(args)
    img_data = all_images.filter { stem, data_map ->
        data_map.data_type in ["raw_image_data", "label_image_data"] // only data to be used to generate images
    }
    .map { stem, data_map -> // structure as Generate_image input
        [
            stem,
            data_map.prefix,
            data_map.data_type.replace("_image_data",""), // make type only `raw` or `label`
            file(data_map.data_path),
            data_map.file_type,
            file(data_map.ref_img ?: "NO_REF") , // set to NO_REF if no ref_img provided. Must be inside a file() to be recognized as a path
            data_map.args ?: [:]
        ]
    }

    Generate_image(img_data)

    Generate_image.out
        .map { stem, prefix, type, paths -> // structure as image_to_zarr input
            [
                stem,
                prefix,
                type,
                [paths].flatten(),
                true // keep_filename, output name was set to stem-prefix-type when generating the image so keep the filename
            ]
        }
        .transpose(by: 3) // make an item for each path of flattened path array like [ stem, prefix, type, path, true ]
        .set {label_tifs}

    all_tifs = img_tifs.mix(label_tifs) // mix img_tifs and label_tifs in no particular order to then put to image_to_zarr
    all_tifs.tap{tifs} // unnecessary (could just use all_tifs)? set assign all_tifs to tifs, used for workflow output
    image_to_zarr(all_tifs)

    ome_zarr_metadata(image_to_zarr.out.ome_xml) // for each image zarr get metadata

    img_zarrs = image_to_zarr.out.img_zarr
        .join(ome_zarr_metadata.out, by: [0,1]) // join zarrs with respective metadata by [ stem, type ]. assuming single raw and single label per dataset ?
        .map { stem, type, path, md ->
            [
                stem, type, [path: path.name, md: new JsonSlurper().parseText(md.trim())] // use JsonSlurper to parse metadata from ome_zarr_metadata output which is just stdout
            ]
        }
        .groupTuple(by: [0,1]) // group by stem and type like [ stem, type, [[path, md], [path, md]] ]. but always only one [path, md] ?
    // provide zarrs and tifs
    emit:
    img_zarrs = img_zarrs // for vitessce config
    img_tifs = tifs // for spatialdata
}


workflow Output_to_config {
    take: 
    out_file_paths
    out_img_zarrs
    
    main:

        // Map workflows' outputs to:
        // tuple val(stem), val(files), val(img_map), val(config_map)

        out_img_zarrs
            .map { stem, type, img -> 
                [stem, [type: type, img: img]] // reduce type and img into a single map with keys `type` and `img`. img is [[path, md]]
            }
            .branch { stem, data -> // make img_zarrs branches for each image type
                raw: data.type == "raw"
                label: data.type == "label"
            }
            .set{img_zarrs}

        img_zarrs.raw
            .join(img_zarrs.label, remainder: true) // join raw images with their respective label (by stem). set remainder: true to allow null values
            .map { stem, raw_data, label_data -> [
                stem,
                [
                    raw: raw_data ? raw_data.img : [], // set raw img or empty
                    label: label_data ? label_data.img : [] // set label img or empty
                ]
            ]}
            .set{img_map}

        datasets.config_map // vitessce config params
            .join(out_file_paths, remainder: true) // join (by stem) with file paths, anndata-zarr and other. allow null
            .join(img_map, remainder: true) // join (by stem) with images. allow null
            .set{data_for_config} // assign to data_for_config
        // items would be like [stem, config_map, files, img_map]

        Build_config(
            data_for_config
            )
}


workflow Output_to_spatialdata {
    take: 
    anndata_files
    img_tifs
    
    main:

        // Similar to Output_to_config
        // Joins (by stem) anndata-zarrs and images into single item in channel to then put through write_spatialdata
        
        img_tifs
            .map { stem, prefix, type, img, k -> // k is keep_filename, ignored in output map
                [stem, [type: type, img: img]]
            }
            .branch { stem, data ->
                raw: data.type == "raw"
                label: data.type == "label"
            }
        .set{tif_files}

        anndata_files
            .join(tif_files.raw, remainder: true)
            .join(tif_files.label, remainder: true)
            .map { stem, anndata, raw_tif, label_tif -> [
                stem, anndata,
                raw_tif ? raw_tif.img : [],
                label_tif ? label_tif.img : []
            ]}
            .set{data_for_sd}

        write_spatialdata(
            data_for_sd
        )
        
}
