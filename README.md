[![python-tests](https://github.com/haniffalab/vitessce-pipeline/actions/workflows/tests-python.yml/badge.svg)](https://github.com/haniffalab/vitessce-pipeline/actions/workflows/tests-python.yml)
[![codecov](https://codecov.io/gh/haniffalab/vitessce-pipeline/branch/main/graph/badge.svg?token=7HQVFH08WJ)](https://codecov.io/gh/haniffalab/vitessce-pipeline/branch/main)

# Vitessce Pipeline

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://haniffalab.github.io/vitessce-pipeline)
[![demo](https://img.shields.io/badge/Demos-view-blue)](https://haniffalab.github.io/vitessce-pipeline/demos.html)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.7405818.svg)](https://doi.org/10.5281/zenodo.7405818)

This Nextflow pipeline processes spatial and single-cell experiment data for visualisation in [vitessce-app](https://github.com/haniffalab/vitessce-app). The pipeline generates data files for [supported data types](http://vitessce.io/docs/data-types-file-types/), and builds a [view config](http://vitessce.io/docs/view-config-json/).


## Usage

The pipeline can handle data from `h5ad` files, SpaceRanger output, Xenium output and MERSCOPE output.

Running the pipeline requires two parameters files that define the data to be processed.
Full instructions and parameters definitions for this files are available in the [documentation](https://haniffalab.com/vitessce-pipeline/setup.html)

One is a `run parameters` yaml file for general configuration options, which can look like

```yaml
data_params: "/path/to/data-parameters.csv"
outdir: "/path/to/output/"
vitessce_options:
  spatial:
    xy: "obsm/spatial"
  mappings:
    obsm/X_umap: [0,1]
  matrix: "X"
layout: "simple"
```

that points to a `data parameters` csv/tsv file which defines files and dataset information, like

```csv
project,dataset,data_type,data_path,args,prefix
project_1,dataset_1,title,Sample Visium Dataset,,
project_1,dataset_1,h5ad,/path/to/visium/anndata.h5ad,,
project_1,dataset_1,raw_image,/path/to/visium/raw_image.tif,,
project_1,dataset_1,label_image,/path/to/visium/label_image.tif,,
project_1,dataset_2,title,Sample Xenium Dataset,,
project_1,dataset_2,xenium,/path/to/xenium/output/,,
project_1,dataset_2,raw_image,/path/to/xenium/raw_image.tif,,
project_1,dataset_2,label_image_data,/path/to/xenium/output/,'{"file_type": "xenium", "shape": [1000,1000]}',
```

The pipeline can then be run like

```sh
nextflow run main.nf -params-file /path/to/run-params.yaml -entry Full_pipeline
```
