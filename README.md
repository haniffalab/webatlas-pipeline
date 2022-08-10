## SINGLE CELL INSIGHTS

# Spatial Data Pipeline

Nextflow pipeline to pre-process spatial data (In-Situ Sequencing, 10x Visium) for [Vitessce](http://github.com/hms-dbmi/vitessce/#readme). The pipeline generates data files for [supported data types](http://vitessce.io/docs/data-types-file-types/), and builds a [view config](http://vitessce.io/docs/view-config-json/).

## Installation

1. Clone the repository

```sh
git clone git@github.com:haniffalab/sci-spatial-data.git
cd sci-spatial-data
```

2. Install nextflow by following the [official instruction](https://www.nextflow.io/index.html#GetStarted)
3. Install [Docker](https://docs.docker.com/engine/install/) and make sure it's in PATH


## Setup

General configuration options are specified with a [yaml file](templates/visium_template.yaml). While dataset-specific configurations are written in a [tsv file](templates/visium_template.tsv)

### yaml file

`outdir` is the path to the directory to which files will be written

`tsv` is the path to the file that contains a dataset to process per line

`args` is a map of arguments for supported files that will be processed. This is applied to files of all datasets. For example,
```yaml
args:
    h5ad:
        compute_embeddings: "True"
    molecules:
        delimiter: "','" 
```

`options` is a map of the contents of the h5ad file that gets converted to Zarr. This is information required to write the Vitessce config file. For example,
```yaml
options:
    spatial:
        xy: "obsm/spatial"
    mappings:
        obsm/X_umap: [0,1]
        obsm/spatial: [0,1]
    factors:
        - "obs/sample"
    sets:
        - "obs/sample"
    matrix: "X"
```
Where `spatial` specifies where the Anndata object holds spatial coordinates. `mappings` is a list of embeddings and the index of the dimensions to be used in a scatterplot. `factors` is a list of useful metadata to be shown per cell when hovering over them in the visualization. `sets` is a list of metadata that defines different cell sets. `matrix` is the expression matrix to use, which will typically be `X`. These values are applied to all datasets unless a dataset has its own `options` specified within the `tsv` file.
**Note** that the pipeline does not check for the existence of these metadata within the h5ad file. It is written directly to the Vitessce config file.

`layout` is the predefined Vitessce layout to use, it can be either `minimal`, `simple` or `advanced`

`custom_layout` is an optional string that defines a Vitessce layout following [Vitessce's View Config API's layout alternative syntax](https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout) (Vitessce components are concatenated horizontally with `|` and vertically with `/`). Supercedes `layout`


### tsv file

Each dataset to be processed is defined as a line in a [tsv file](templates/visium_template.tsv)

`h5ad` is the path to the h5ad file

`image_type` specifies whether the image at `image_path` is type `raw` or `label`

`title` is the project/experiment title which can have multiple datasets

`dataset` is the name of the dataset

`molecule` is the path to a molecules tsv file

`url` is an optional string that will be prepended to each file in the Vitessce config file

`options` is an optional json-like string that overrides the `options` in the yaml file


### Docker

Before running the pipeline, build the docker images.

```sh
cd docker
./build-docker-imgs.sh
```

## Run

The pipeline contains workflows to process files (h5ad files, csv/tsv files), convert images to Zarrs and build a Vitessce config file from the generated files.

The `scRNAseq_pipeline` entry point handles datasets with no image data. This processes files and builds a Vitessce config file.
The `ISS_pipeline` entry point handles datasets that contain both raw and label images.
The `Visium_pipeline` entry point handles datasets that contain a raw image and generates a label image from the specified h5ad file where it expects to have a `spatial` key within `uns`.

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

Further reading:
--- 

Docker image pulling/local conda env creation are handled by nextflow. Please refer to [this](https://www.nextflow.io/docs/latest/getstarted.html) for detailed information.


## Testing

### Python testing

Testing of python scripts uses [pytest](https://docs.pytest.org/en/7.1.x/).

Set the `PYTHONPATH` environment variable to the `bin` directory where the scripts are stored and then run
```
python -m pytest -q tests/test_class.py
```
