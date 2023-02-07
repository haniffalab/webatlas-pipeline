.. _setup:

#####
Setup
#####

Currently, the pipeline can process several types of data files as well as images.
It can handle ``h5ad`` files, SpaceRanger output, Xenium output, MERSCOPE output, and ``molecules`` data in ``csv``/``tsv`` files.
It can handle raw and label images in ``tif`` format, as well as generating/preprocessing images from Visium, Xenium and MERFISH data. 

Running the pipeline requires two parameters files:

1. a `run parameters`_ ``yaml`` file with general configuration options
2. a `data parameters`_ ``csv``/``tsv`` file with dataset-specific information

Templates are available in the `templates directory <templates/>`__.

.. _run-parameters:

**************
Run parameters
**************

The run parameters file is a ``yaml`` that defines several configuration options which apply to a full run of the pipeline (which can process one or more datasets).
Parameters for each type of data file work as the default parameters to all the data files defined in the `data parameters`_ file.

In the command line, the path to the ``yaml`` file is indicated through the ``-params-file`` flag.

A run parameters file can look like:

.. code-block:: yaml

    data_params: "/path/to/data-parameters.csv"
    outdir: "/path/to/output/"
    args:
        h5ad:
            compute_embeddings: "True"
    vitessce_options:
        spatial:
            xy: "obsm/spatial"
        mappings:
            obsm/X_umap: [0,1]
    layout: "simple"


The parameters are as follows:

- ``data_params``, the path to the `data parameters`_ ``csv``/``tsv`` file that contains the datasets' information and data files to process.

- ``data_params_delimiter``, the delimiter to be used for the `data parameters`_ file.
  Defaults to ``,``.

- ``outdir``, the path to the directory to which output files will be written.

- ``args``, a map of optional arguments for the scripts that process the supported files. 
  This is applied to files of all datasets. See available `args`_.

- ``vitessce_options``, a map of the contents of a dataset's Anndata object that gets
  converted to Zarr. This is information required to write the Vitessce
  config file. These values are applied to all datasets unless a dataset
  has its own ``options`` specified within the `data parameters`_ file. 
  :ref:`See more <run-parameters-vitessce-options>`.

- ``layout``, a predefined Vitessce layout to use, it can be either
  ``minimal``, ``simple`` or ``advanced``.

- ``custom_layout``, an optional string that defines a Vitessce layout
  following `Vitessce's View Config API's layout alternative
  syntax <https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout>`__
  (Vitessce components are concatenated horizontally with ``|`` and
  vertically with ``/``). Supersedes ``layout``.

.. _run-parameters-args:

Args
====

Available ``args`` depend of the supported data type of each data file.

Supported data types are ``h5ad``, ``spaceranger``, ``xenium``, ``merscope``, ``molecules``

This are passed to the processing script for all the files of the 
corresponding data type (unless overriden in the `data parameters`_ file).
Possible values for the currently supported data types are as follows:

.. code-block:: yaml

  args:
    h5ad:
      compute_embeddings: "True" # set to `True` to compute PCA and UMAP if not already within the anndata object
      chunk_size: 20 # Zarr chunk size, defaults to 10
      var_index: "SYMBOL" # `var` column from the anndata object to use as the gene names in the webapp. This reindexes the `var` matrix
      obs_subset: ["sample", ["sample_id_1"]] # optional `obs` column name an value(s) to subset the anndata object
      var_subset: ["genome", ["GRCh38"]] # optional `var` column name an value(s) to subset the anndata object
      batch_processing: "False" # set to `True` to process the file in batches to avoid loading the whole object into memory if it is too large
      batch_size: 1000 # batch size (number of columns to process at a time if matrix is dense/csc, number of rows if matrix is csr) if `batch_processing` is set to `True`
    spaceranger:
      save_h5ad: "True" # save the intermediate h5ad to the output directory. Defaults to `False`
      load_clusters: "True" # set to `False` to disable loading the clusters from the `analysis` directory
      load_embeddings: "True" # set to `False` to disable loading the embeddings (UMAP, tSNE and PCA) from the `analysis` directory
      clustering: "graphclust" # clusters to load (if `load_clusters` is set to `True`)
    xenium:
      save_h5ad: "True" # save the intermediate h5ad to the output directory. Defaults to `False`
      spatial_as_pixel: "True" # convert spatial coordinates to pixel coordinates. Defaults to `True`
      resolution: 0.2125 # pixel resolution used to convert the spatial coordinates. Defaults to 0.2125
    merscope:
      save_h5ad: "True" # save the intermediate h5ad to the output directory. Defaults to `False`
      filter_prefix: "Blank-" # prefix to filter out data from the cell by gene data. Defaults to `Blank-`
    molecules:
      delimiter: "," # the file delimiter. Defaults to `\t`
      has_header: "True" # set to `False` if csv/tsv file contains no header
      gene_col_name: "Name" # name of the column for gene names. Defaults to `Name`.
      x_col_name: "x_int" # name of the column for `x` coordinates. Defaults to `x_int`.
      y_col_name: "y_int" # name of the column for `y` coordinates. Defaults to `y_int`.
      gene_col_idx: 0 # column index of the column for gene names in case `has_header` is `False`.
      x_col_idx: 1 # column index of the column for `x` coordinates in case `has_header` is `False`.
      y_col_idx: 2 # column index of the column for `y` coordinates in case `has_header` is `False`.

Note that in the case of ``spaceranger``, ``xenium`` and ``merscope`` data, it initially gets converted 
into an ``h5ad`` file and so when processed the ``args`` for ``h5ad`` also apply to it, 
unless overriden in the `data parameters`_ file.

.. _run-parameters-vitessce-options:

Vitessce options
================

The ``vitessce_options`` map is used to write Vitessce config files.
One Vitessce config file is generated per dataset.
Include relevant information from your data to be visualized.
All values are optional as they depend on them existing in your data.

Values that can be specified are ``spatial``, ``mappings``, ``factors``, ``sets`` and ``matrix``,
and should be defined as the following example:

.. code-block:: yaml

  vitessce_options:
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

**Note** that the pipeline does not check for the existence of these
metadata within the h5ad file. It is written directly to the Vitessce
config file. If they're incorrectly specified then an error will occur when
Vitessce tries to load the data. The output config
file can be manually edited without re-running the pipeline to fix or adapt 
the visualization to your needs.

.. _data-parameters:

***************
Data parameters
***************

The data parameters file is a ``csv``/``tsv`` file used to define `dataset information`_ and :ref:`data files <data-parameters-data-file>` to be processed.
Multiple datasets can be defined in the same data parameters file and they will all be processed in the same pipeline run.
Each line can either define a file/image *or* dataset information.

For example:

.. code-block:: text

    project,dataset,data_type,data_path,args,prefix
    project_1,dataset_1,title,Sample dataset,,
    project_1,dataset_1,h5ad,/path/to/visium/anndata.h5ad,,
    project_1,dataset_1,raw_image,/path/to/raw_image.tif,,
    project_1,dataset_1,label_image,/path/to/label_image.tif,,


The supported image format is ``tif``.
Images can be either raw images (microscopy images) or label images (containing segmentations).
Additionally, label images can be generated and processed if provided with the necessary data.
Label images can be generated for Visium data if provided with an ``h5ad`` file or
SpaceRanger output directory, and Xenium and MERSCOPE if provided with their 
respective output directories.
Raw images can also be pre-processed, in the case of MERSCOPE data where the raw image channels
are stored in separate ``tif`` files the pipeline can concatenate them to then convert them.

Datasets are identified and grouped by the joint ``project-dataset`` key.
Each dataset does not need to contain all types of data types,
but it should contain at least one file or image be processed.

Columns definitions:

- ``project`` is the project/experiment name which can have multiple datasets

- ``dataset`` is the name of the dataset

- ``data_type`` is the type of file to be processed if the line defines a `data file`_,
  or the type of `dataset information`_ otherwise.

- ``data_path`` is the path to file or directory containing the data if the line defines a `data file`_,
  or the `dataset information`_ value.

- ``args`` is an optional JSON-like string defining argument names and values 
  to be used in the script that processes the data file.
  It must be written inside simple quotes ``'``, with strings inside it using double quotes ``"``,
  like ``'{"key": "value"}'``.
  This overrides ``args`` from `run parameters`_ for the line's file only.
  This value is not used if the line is defining dataset information.

- ``prefix`` is an optional string to prefix the output filenames, along with the ``project``
  and ``dataset`` names, so the output filenames become ``{project}-{dataset}-{prefix}-file.ext``.
  Required if you have multiple input files of the same ``data_type`` within the same ``project``
  and ``dataset``, as they would otherwise get 
  overwritten with the default output filename ``{project}-{dataset}-file.ext``.
  If a single input file generates multiple output files of the same type, a prefix will
  automatically be added to each of them to avoid overwritting.

.. _data-parameters-data-file:

Data file
=========

A line defining a data file can be written as follows:


.. code-block:: text

    project,dataset,data_type,data_path,args,prefix
    project_1,dataset_1,h5ad,/path/to/visium/anndata.h5ad,,

Supported values are 

.. list-table:: 
    :widths: 10 10 15
    :header-rows: 1

    * - data_type
      - data_path
      - args
    * - ``h5ad``
      - Path to the ``h5ad`` file
      - JSON-like string with arguments as described in `args`_
    * - ``spaceranger``
      - Path to a SpaceRanger output directory
      - JSON-like string with arguments as described in `args`_
    * - ``xenium``
      - Path to a Xenium output directory
      - JSON-like string with arguments as described in `args`_
    * - ``merscope``
      - Path to a MERSCOPE output directory
      - JSON-like string with arguments as described in `args`_
    * - ``molecules``
      - Path to a molecules ``csv``/``tsv`` file
      - JSON-like string with arguments as described in `args`_
    * - ``raw_image``
      - Path to the raw ``tif`` image
      - None
    * - ``label_image``
      - Path to the raw ``tif`` image
      - None
    * - ``raw_image_data``
      - Path to a file or directory containing data from which to generate or pre-process a raw ``tif`` image.
        
        Possible inputs depend on the supported technology from which the data is obtained,
          
          * ``merscope`` requires a path to the output directory containing an ``images`` directory
            where image channels are stored as ``tif`` files
      - JSON-like string with data type-specific key-values. 
        
        :ref:`See full list <data-parameters-raw-image-data>`
    * - ``label_image_data``
      - Path to a file or directory containing data from which to generate a label ``tif`` image. 
        
        Possible inputs depend on the supported technology from which the data is obtained,
          
          * ``visium`` requires a path to an ``h5ad`` file or SpaceRanger output directory
          * ``xenium`` requires a path to a Xenium output directory
          * ``merscope`` requires a path to a MERSCOPE output directory
      - JSON-like string with data type-specific key-values.
        
        :ref:`See full list <data-parameters-label-image-data>`


.. _data-parameters-raw-image-data:

raw image data
^^^^^^^^^^^^^^

Supported raw image data types: ``merscope``

``args`` are defined with a JSON-like string containing the key-values corresponding to the type of data provided
as described below

.. list-table:: 
    :widths: 10 10 15
    :header-rows: 1

    * - type
      - key
      - value
    * - all
      - ``file_type`` 
        (required)
      - ``merscope``
    * - ``merscope``
      - ``z_index`` 
        (optional)
      - Z index or list of Z indices of which to concatenate their separate channel images into a single ``tif`` image

        Defaults to ``[0]``

For example:

.. code-block:: text

    project,dataset,data_type,data_path,args,prefix
    project_1,dataset_1,raw_image_data,/path/to/merscope/output/,'{"file_type": "merscope"}',
    project_1,dataset_2,raw_image_data,/path/to/merscope2/output/,'{"file_type": "merscope", "z_index": [0,1]}',


.. _data-parameters-label-image-data:

label image data
^^^^^^^^^^^^^^^^

Supported label image data types: ``visium``, ``xenium``, ``merscope``

``args`` are defined with a JSON-like string containing the key-values corresponding to the type of data provided
as described below

.. list-table::
    :widths: 10 10 15
    :header-rows: 1

    * - type
      - key 
      - value
    * - all
      - ``file_type``
        (required)
      - ``visium``, ``xenium`` or ``merscope``.
    * - all
      - ``ref_img``
        (optional)
      - A reference ``tif`` image of the size of the desired label image.
    * - all
      - ``shape``
        (optional)
      - Shape of the desired label image as ``[int, int]``.
    * - ``visium``
      - ``obs_subset``
        (optional)
      - A tuple containing ``obs`` column name and value(s) like to subset the Anndata object.
        Useful if the the Anndata object contains data of multiple slides.
    * - ``visium``
      - ``sample_id``
        (optional)
      - The name of the sample within the Anndata object.
        Otherwise the first one will be used.
    * - ``xenium``
      - ``resolution``
        (optional)
      - Pixel resolution.
        
        Defaults to ``0.2125``
    * - ``merscope``
      - ``z_index``
        (optional)
      - Z index or list of Z indices of label images to generate.

        Defaults to ``[0]``

For example:

.. code-block:: text

    project,dataset,data_type,data_path,args,prefix
    project_1,dataset_1,label_image_data,/path/to/visium/anndata.h5ad,'{"file_type": "visium", "ref_img": "/path/to/raw_image.tif", "sample_id": "sample_1"}',
    project_1,dataset_2,label_image_data,/path/to/xenium/output/,'{"file_type": "xenium", "shape": [1000,1000]}',
    project_1,dataset_3,label_image_data,/path/to/merscope/output/,'{"file_type": "merscope", "ref_img": "/path/to/raw_image2.tif"}',

.. _data-parameters-dataset-info:

Dataset information
===================

A line defining optional dataset information can be written as follows

.. code-block:: text

    project,dataset,data_type,data_path,args,prefix
    project_1,dataset_1,title,Dataset 1,,

Supported values are 

.. list-table:: 
    :widths: 10 15
    :header-rows: 1

    * - data_type
      - data_path
    * - ``title``
      - Name or title for the final Vitessce config file and visualization.
    * - ``description``
      - Dataset description 
    * - ``url``
      - The url to prepend to each converted data file in the output Vitessce config file.
        Vitessce will load files from this location.
        This may be the final location to which files will be uploaded to and served
        or a local one for testing. 
        
        Defaults to ``http://localhost:3000/``
    * - ``layout``
      - A predefined Vitessce layout to use, it can be either 
        ``minimal``, ``simple`` or ``advanced``.
        
        Overrides ``layout`` from `run parameters`_.
    * - ``custom_layout``
      - An optional string that defines a Vitessce layout
        following `Vitessce's View Config API's layout alternative
        syntax <https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout>`__.
        
        Overrides ``custom_layout`` from `run parameters`_.
    * - ``vitessce_options``
      - (*Not recommended*) JSON-like string of values as described in `vitessce options`_.
        This will override the ``vitessce_options`` defined in `run parameters`_ for a specific
        dataset only. Though, the numerous values needed would result in a lengthy string,
        therefore we **strongly recommend** writing another `run parameters`_ file instead of 
        overriding ``vitessce_options``.

Note no ``args`` or ``prefix`` are required for any type of dataset information.

.. _setup-docker :

******
Docker
******

Before running the pipeline, build the docker images.

.. code-block:: sh

   cd docker
   ./build-docker-imgs.sh
