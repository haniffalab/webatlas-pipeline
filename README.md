## SINGLE CELL INSIGHTS

[![tests](https://github.com/haniffalab/sci-spatial-data/actions/workflows/tests-python.yml/badge.svg)](https://github.com/haniffalab/sci-spatial-data/actions/workflows/tests-python.yml)
[![codecov](https://codecov.io/gh/haniffalab/sci-spatial-data/branch/dev/graph/badge.svg?token=7HQVFH08WJ)](https://codecov.io/gh/haniffalab/sci-spatial-data/branch/dev)
[![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)

# Spatial Data Pipeline

Nextflow pipeline to pre-process spatial data (In-Situ Sequencing, 10x Visium) for [Vitessce](http://github.com/hms-dbmi/vitessce/#readme). (Single cell RNA sequencing data can also be visualized without the spatial layer.) The pipeline generates data files for [supported data types](http://vitessce.io/docs/data-types-file-types/), and builds a [view config](http://vitessce.io/docs/view-config-json/).

# Installation

1. Clone the repository

```sh
git clone git@github.com:haniffalab/sci-spatial-data.git
cd sci-spatial-data
```

2. Install nextflow by following the [official instruction](https://www.nextflow.io/index.html#GetStarted)
3. Install [Docker](https://docs.docker.com/engine/install/) and make sure it's in PATH


# Setup

General configuration options are specified with a [yaml file](templates/visium_template.yaml). While dataset-specific configurations are written in a [tsv file](templates/visium_template.tsv)

Templates are available in the [templates directory](templates/).

## yaml file

`outdir` is the path to the directory to which files will be written

`tsv` is the path to the file that contains a dataset to process per line

`args` is a map of optional arguments for the scripts that process the supported files. This is applied to files of all datasets.
Available `args` are
```yaml
args:
    h5ad:
        compute_embeddings: "True" # set to `True` to compute PCA and UMAP if not already within the anndata object
        chunk_size: 20 # Zarr chunk size, defaults to 10
        var_index: "SYMBOL" # `var` column from the anndata object to use as the gene names in the webapp. This reindexes the `var` matrix
    spaceranger:
        save_h5ad: "True" # save an h5ad to the output directory. Defaults to `False`
        load_clusters: "True" # set to `False` to disable loading the clusters from the `analysis` directory
        load_embeddings: "True" # set to `False` to disable loading the embeddings (UMAP, tSNE and PCA) from the `analysis` directory
        clustering: "graphclust" # clusters to load (if `load_clusters` is set to `True`)
    molecules:
        delimiter: "','" # the file delimiter. Defaults to `\t`
        has_header: "True" # set to `False` if csv/tsv file contains no header
        gene_col_name: "Name" # name of the column for gene names. Defaults to `Name`.
        x_col_name: "x_int" # name of the column for `x` coordinates. Defaults to `x_int`.
        y_col_name: "y_int" # name of the column for `y` coordinates. Defaults to `y_int`.
        gene_col_idx: 0 # column index of the column for gene names in case `has_header` is `False`.
        x_col_idx: 1 # column index of the column for `x` coordinates in case `has_header` is `False`.
        y_col_idx: 2 # column index of the column for `y` coordinates in case `has_header` is `False`.
```

`options` is a map of the contents of the h5ad file that gets converted to Zarr. This is information required to write the Vitessce config file. These values are applied to all datasets unless a dataset has its own `options` specified within the `tsv` file. For example,
```yaml
options:
    spatial:
        xy: "obsm/spatial" # where the Anndata object holds spatial coordinates
    mappings: # list of embeddings and the index of the dimensions to use in a scatterplot
        obsm/X_umap: [0,1]
        obsm/spatial: [0,1]
    factors: # list of useful metadata to show per cell when hovering over them in the visualization
        - "obs/sample"
    sets: # list of keys for grouping cells
        - name: "obs/celltype" # key with cell set labels
          score: "obs/celltype_prob" # key with cell set confidence/percentage scores (float values [0,1])
        - "obs/sample" # key with cell set labels, without associated scores
    matrix: "X" # expression matrix to use
```
**Note** that the pipeline does not check for the existence of these metadata within the h5ad file. It is written directly to the Vitessce config file. If they're incorrectly specified then the output config file can be manually edited without re-running the pipeline.

`layout` is the predefined Vitessce layout to use, it can be either `minimal`, `simple` or `advanced`

`custom_layout` is an optional string that defines a Vitessce layout following [Vitessce's View Config API's layout alternative syntax](https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout) (Vitessce components are concatenated horizontally with `|` and vertically with `/`). Supercedes `layout`


## tsv file

Each dataset to be processed is defined as a line in a [tsv file](templates/visium_template.tsv)

`h5ad` is the path to the h5ad file

`spaceranger` is the path to the spaceranger output directory

`image_type` specifies whether the image at `image_path` is type `raw` or `label`

`title` is the project/experiment title which can have multiple datasets

`dataset` is the name of the dataset

`molecule` is the path to a molecules tsv file

`url` is an optional string that will be prepended to each file in the Vitessce config file

`args` is an optional json-like string that overrides the `args` in the yaml file

`options` is an optional json-like string that overrides the `options` in the yaml file


## Docker

Before running the pipeline, build the docker images.

```sh
cd docker
./build-docker-imgs.sh
```

# Run

The pipeline contains workflows to process files (h5ad files, spaceranger output, csv/tsv files), convert images to Zarrs and build a Vitessce config file from the generated files.

- The `scRNAseq_pipeline` entry point handles datasets with no image data. This processes files and builds a Vitessce config file.
- The `ISS_pipeline` entry point handles datasets that contain both raw and label images.
- The `Visium_pipeline` entry point handles datasets that contain a raw image and generates a label image from the specified h5ad file where it expects to have a `spatial` key within `uns`.

Each dataset in a tsv file should belong to the same modality so they can all be run through the corresponding pipeline. 

```
nextflow run main.nf -params-file [your_params].yaml -entry [modality]_pipeline
```

You can also just convert the images to Zarrs:

```
nextflow run main.nf -params-file [your_params].yaml -entry To_ZARR
```

Or just process files:

```
nextflow run main.nf -params-file [your_params].yaml -entry Process_files
```

Additionally, you can build a config file without processing files nor images, using a slightly different [yaml file](templates/config_template.yaml), per dataset, where you would provide filenames and necessary metadata:

```
nextflow run main.nf -params-file [your_params].yaml -entry Config
```

## Further reading:

Docker image pulling/local conda env creation are handled by nextflow. Please refer to [this](https://www.nextflow.io/docs/latest/getstarted.html) for detailed information.


# Testing

## Python testing

Testing of python scripts uses [pytest](https://docs.pytest.org/en/7.1.x/).

Set the `PYTHONPATH` environment variable to the `bin` directory where the scripts are stored and then run
```
python -m pytest -q tests/test_class.py
```
