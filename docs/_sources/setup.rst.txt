.. _setup:

Setup
=====

Running the pipeline requires two parameters files:

1. a `run parameters`_ ``yaml`` file with general configuration options
2. a `data parameters`_ ``csv``/``tsv`` file with dataset-specific information

Templates are available in the `templates directory <templates/>`__.

.. _run-parameters:

Run parameters
--------------

The parameters in the ``yaml`` file apply to a full run of the pipeline.
Parameters for each type of data file work as the default parameters to all the data files defined in the `data parameters`_ file.

In the command line, the path to the ``yaml`` file is indicated through the ``-params-file`` flag.

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
  See more information about `vitessce options`_.

- ``layout``, a predefined Vitessce layout to use, it can be either
  ``minimal``, ``simple`` or ``advanced``.

- ``custom_layout``, an optional string that defines a Vitessce layout
  following `Vitessce's View Config API's layout alternative
  syntax <https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout>`__
  (Vitessce components are concatenated horizontally with ``|`` and
  vertically with ``/``). Supersedes ``layout``.

.. _run-parameters-args:

args
^^^^

Available ``args`` depend of the data type of each data file.
This are passed to the processing script for all the files of the 
corresponding data type (unless overriden in the `data parameters`_ file).
Possible values for the currently supported data types are as follows:

.. code:: yaml

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
      delimiter: "," # the file delimiter. Defaults to `\t`
      has_header: "True" # set to `False` if csv/tsv file contains no header
      gene_col_name: "Name" # name of the column for gene names. Defaults to `Name`.
      x_col_name: "x_int" # name of the column for `x` coordinates. Defaults to `x_int`.
      y_col_name: "y_int" # name of the column for `y` coordinates. Defaults to `y_int`.
      gene_col_idx: 0 # column index of the column for gene names in case `has_header` is `False`.
      x_col_idx: 1 # column index of the column for `x` coordinates in case `has_header` is `False`.
      y_col_idx: 2 # column index of the column for `y` coordinates in case `has_header` is `False`.

.. _run-parameters-vitessce-options:

vitessce options
^^^^^^^^^^^^^^^^

The ``vitessce_options`` map is used to write the Vitessce config file.
One Vitessce config file is generated per dataset.
Include relevant information from your data to be visualized.
All values are optional as they depend on them existing in your data.

Values that can be specified are ``spatial``, ``mappings``, ``factors``, ``sets`` and ``matrix``,
and should be defined as the following example:

.. code:: yaml

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

Data parameters
---------------

The ``csv``/``tsv`` file is used to define dataset information and data files to be processed.
Multiple datasets can be defined in the same `data parameters` file and they will all be processed in the same pipeline run.
Datasets are identified and grouped by the joint ``project-dataset`` key.

Each dataset does not need to contain all types of data types,
but it should contain at least one to be processed (file or image).


Each line can either define a data file *or* dataset information.
Examples for each case are provided further down.

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

.. _data-parameters-data-file:

data file
^^^^^^^^^

A line defining a data file can be written as follows::

    project,dataset,data_type,data_path,args
    project_1,dataset_1,h5ad,/path/to/visium/anndata.h5ad,

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
    * - ``molecules``
      - Path to a molecules ``csv``/``tsv`` file
      - JSON-like string with arguments as described in `args`_
    * - ``raw_image``
      - Path to the raw ``tif`` image
      - None
    * - ``label_image``
      - Path to the raw ``tif`` image
      - None
    * - ``label_image_data``
      - Path to a file or directory containing data from which to generate a label ``tif`` image. 
        
        Possible inputs depend on the supported technology from which the data is obtained,
          * ``visium`` requires a path to an ``h5ad`` file or ``spaceranger`` output directory
      - JSON-like string with the following key-values,
          * ``file_type`` (required), supported technology like ``visium``.
          * ``ref_img`` (optional), a reference ``tif`` image of the size of the desired label image
          * ``shape`` (optional), shape of the desired label image as ``[int, int]``
          * ``sample_id`` (optional), the name of the sample within the Anndata object.
            Otherwise the first one will be used.

        For example,

        ``'{"file_type": "visium", "ref_img": "/path/to/raw.tif", "sample_id": "visium_sample"}'``

        or

        ``'{"file_type": "visium", "shape": [1000,1000], "sample_id": "visium_sample"}'``

.. _data-parameters-dataset-info:

dataset information
^^^^^^^^^^^^^^^^^^^

A line defining optional dataset information can be written as follows::

    project,dataset,data_type,data_path,args
    project_1,dataset_1,url,http://localhost:3000/visium_dataset_1/,

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
    * - ``layout``
      - a predefined Vitessce layout to use, it can be either 
        ``minimal``, ``simple`` or ``advanced``.
        Overrides ``layout`` from `run parameters`_.
    * - ``custom_layout``
      - an optional string that defines a Vitessce layout
        following `Vitessce's View Config API's layout alternative
        syntax <https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout>`__.
        Overrides ``custom_layout`` from `run parameters`_.
    * - ``vitessce_options``
      - (*Not recommended*) JSON-like string of values as described in `vitessce options`_.
        This will override the ``vitessce_options`` defined in `run parameters`_ for a specific
        dataset only. Though, the numerous values needed would result in a lengthy string,
        therefore we **strongly recommend** writing another `run parameters`_ file instead of 
        overriding ``vitessce_options``.

Note no ``args`` are required for any type of dataset information.

.. _setup-docker :

Docker
------

Before running the pipeline, build the docker images.

.. code:: sh

   cd docker
   ./build-docker-imgs.sh
