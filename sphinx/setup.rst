.. _setup:

#####
Setup
#####

Currently, the pipeline can process several types of data files as well as images.

- It can handle ``h5ad`` files, SpaceRanger output, Xenium output, MERSCOPE output, and ``molecules`` data in ``csv``/``tsv`` files.
- It can handle raw and label images in ``tif`` format, as well as generating/preprocessing images from Visium, Xenium and MERFISH data. 

Running the pipeline requires a ``yaml`` `parameters file`_ that lists the data to be processed.

Templates of this parameters file are available in the `templates directory <templates/>`__.


.. _parameters_file:

***************
Parameters file
***************

The parameters file defines several configuration options which apply to a full run of the pipeline (which can process one or more datasets).

Running the pipeline, the path to the ``yaml`` file is indicated through the ``-params-file`` flag in the command line.

A parameters file looks like:

.. code-block:: yaml

    outdir: "/path/to/output/"
    
    args:
        h5ad:
            compute_embeddings: "True"
    
    projects:
      - project: project_1
        datasets:
          - dataset: dataset_1
            data:
              -
                data_type: h5ad
                data_path: /path/to/project_1/dataset_1/anndata.h5ad

    vitessce_options:
        spatial:
            xy: "obsm/spatial"
        mappings:
            obsm/X_umap: [0,1]
    layout: "simple"
    # custom_layout: "spatial|scatterplot" # overrides `layout`


The parameters are as follows:

.. list-table:: 
    :widths: 10 15
    :header-rows: 1

    * - key
      - value 
    * - ``outdir``
      - the path to the directory to which output files will be written.
    * - ``args``
      - a map of optional arguments for the scripts that process the supported files. 
        
        These serve as the default values for all files in the parameters file. 
        
        See `args`_.
    * - ``projects``
      - a list of projects where each contain a list of datasets that hold one or more
        data files to process. 
        
        See `project`_.
    * - ``vitessce_options``
      - a map of the contents of an input Anndata object
        to be shown in the Vitessce visualization.
        
        These values serve as the default for all projects. 
        
        See `vitessce_options`_.
    * - ``layout``
      - a predefined Vitessce layout to use, it can be either
        ``minimal``, ``simple`` or ``advanced``. 
        
        This value serves as the default for all projects.
    * - ``custom_layout``
      - `optional` string that defines a Vitessce layout
        following `Vitessce's View Config API's layout alternative
        syntax <https://vitessce.github.io/vitessce-python/api_config.html#vitessce.config.VitessceConfig.layout>`__
        (Vitessce components are concatenated horizontally with ``|`` and
        vertically with ``/``). 
        
        Supersedes ``layout``. 
        
        This value serves as the default for all projects.


.. _project:

Project
========

.. code-block:: yaml

  projects:
    - project: project_1
      datasets: 
        ...
      args: 
        h5ad:
          batch_processing: True


Multiple projects can be defined in a single parameters file, and each can define multiple datasets.
Each project item is defined by the following keys

.. list-table:: 
    :widths: 10 15
    :header-rows: 1

    * - key
      - value 
    * - ``project``
      - a unique project name/id
    * - ``datasets``
      - a list of dataset items.
        
        See `dataset`_.
    * - ``args``
      - `optional` map of arguments to set as default for all files within the project.
        
        Supersedes global ``args``. 


.. _dataset:

Dataset
-------

.. code-block:: yaml

  ...
    datasets: 
      - dataset: dataset_1
        title: "Dataset 1"
        data:
          ...
        args:
          h5ad:
            compute_embeddings: True 


Multiple datasets belong to a single `project`_, and each dataset can contain multiple data files.
Each dataset item is defined by the following keys

.. list-table:: 
    :widths: 10 15
    :header-rows: 1

    * - key
      - value 
    * - ``dataset``
      - a unique (within the project) name/id
    * - ``data``
      - a list of data files to be processed. See `data`_.
    * - ``title``
      - `optional` title to show in the Vitessce visualization
    * - ``description`` 
      - `optional` description to show in the Vitessce visualization
    * - ``url``
      - `optional` url to prepend to each converted data file in the output Vitessce config file.
        Vitessce will load files from this location.
        This may be the final location to which files will be uploaded to and served
        or a local location for testing.
        
        Defaults to ``http://localhost:3000/``
    * - ``layout``
      - `optional` predefined Vitessce layout to use.
        
        Supersedes the global ``layout``.
    * - ``custom_layout``
      - `optional` string that defines a Vitessce layout.
        
        Supersedes the global ``custom_layout``. Supersedes ``layout``.
    * - ``vitessce_options``
      - `optional` map of the contents of an input Anndata object
        to be shown in the Vitessce visualization.
        
        Supersedes the global ``vitessce_options``.
    * - ``args``
      - `optional` map of arguments to set as default for all files within the dataset.
        
        Supersedes global and project ``args``.


.. _data:

Data
^^^^

.. code-block:: yaml

  ...
    data:
      - 
        data_type: h5ad
        data_path: /path/to/anndata.h5ad
        args:
          compute_embeddings: True
          batch_processing: True
      - 
        data_type: raw_image
        data_path: /path/to/sample_1.tif
        prefix: sample_1


Each data item defines an input file to be processed by the pipeline.
Multiple data files belong to a single `dataset`_.

A data item can define data files or images.

Data files are processed into AnnData objects and written to Zarr.
Image ilres are written to Zarr through `bioformats2raw`.

For image files the supported image format is ``tif``.
Images can be either raw images (microscopy images) or label images (containing segmentations).
Additionally, label images can be generated and processed if provided with the necessary data.
Label images can be generated for Visium data if provided with an ``h5ad`` file or
SpaceRanger output directory, and Xenium and MERSCOPE if provided with their 
respective output directories.
Raw images can also be pre-processed, in the case of MERSCOPE data where the raw image channels
are stored in separate ``tif`` files the pipeline can concatenate them to then convert them.

Each data item must define at least
its ``data_type`` and ``data_path``.  

Supported values are 

.. list-table:: 
    :widths: 10 10
    :header-rows: 1

    * - data_type
      - data_path
    * - ``h5ad``
      - Path to the ``h5ad`` file
    * - ``spaceranger``
      - Path to a SpaceRanger output directory
    * - ``xenium``
      - Path to a Xenium output directory
    * - ``merscope``
      - Path to a MERSCOPE output directory
    * - ``molecules``
      - Path to a molecules ``csv``/``tsv`` file
    * - ``raw_image``
      - Path to the raw ``tif`` image
    * - ``label_image``
      - Path to the raw ``tif`` image
    * - ``raw_image_data``
      - Path to a file or directory containing data from which to generate or pre-process a raw ``tif`` image.
        
        Possible inputs depend on the supported technology from which the data is obtained,
          
          * ``merscope`` requires a path to the output directory containing an ``images`` directory
            where image channels are stored as ``tif`` files
    * - ``label_image_data``
      - Path to a file or directory containing data from which to generate a label ``tif`` image. 
        
        Possible inputs depend on the supported technology from which the data is obtained,
          
          * ``visium`` requires a path to an ``h5ad`` file or SpaceRanger output directory
          * ``xenium`` requires a path to a Xenium output directory
          * ``merscope`` requires a path to a MERSCOPE output directory


Each data item is defined with the following keys

.. list-table:: 
    :widths: 10 15
    :header-rows: 1

    * - key
      - value 
    * - ``data_type``
      - one of the supported types of files to be processed.
        
        A type of data file: ``h5ad``, ``spaceranger``, ``molecules``, ``xenium``, ``merscope``,
        
        a type of image file ``raw_image``, ``label_image``,
        
        or a type of image/data to generate/preprocess an image ``raw_image_data``, ``label_image_data`` 
    * - ``data_path``
      - path to the file or directory containing the data
    * - ``prefix``
      - `optional` string to prefix the output filenames, along with the ``project``
        and ``dataset`` names, so the output filenames become ``{project}-{dataset}-{prefix}-file.ext``.
        
        Required if you have multiple input files of the same ``data_type`` within the same ``project``
        and ``dataset``, as they would otherwise get 
        overwritten with the default output filename ``{project}-{dataset}-file.ext``.
        
        If a single input file generates multiple output files of the same type, a prefix will
        automatically be added to each of them to avoid overwritting.
    * - ``args``
      - `optional` map of arguments to use when processing the data file.
        
        **Note** that this must be a map of only arguments that correspond to the file's ``data_type``.
        
        Supersedes global, project and dataset ``args`` for that ``data_type``.

In the case where ``data_type`` is ``raw_image_data`` or ``label_image_data``
extra keys should be defined

.. list-table:: 
    :widths: 10 10 15
    :header-rows: 1

    * - key
      - data_type
      - value 
    * - ``file_type``
      - ``raw_image_data`` 
        
        or 
        
        ``label_image_data``
      - ``visium``, ``xenium`` or ``merscope``.
    * - ``ref_img``
      - ``label_image_data``
      - (required if ``shape`` is not set) 
        a reference ``tif`` image of the size of the desired label image.
    * - ``shape``
      - ``label_image_data``
      - (required if ``ref_img`` is not set) 
        shape of the desired label image as ``[int, int]``.


.. _args:

Args
====

Available ``args`` depend of the ``data_type`` of each `data`_ item.

Image files data types ``raw_image`` and ``label_image`` take no ``args``.
Image data files data types (files from which to generate image files or images that require preprocessing)
``raw_image_data`` and ``label_image_data`` take ``args`` depending on their ``file_type`` (``visium``, ``xenium`` or ``merscope``).

Possible values for each of the supported data types are as follows

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
into an ``h5ad`` file and so when processed the ``args`` for ``h5ad`` also apply to them.
If specifying ``args`` directly to a `data`_ item of these types
you can define both the ``args`` for that specific ``data_type`` and ``h5ad``.

Example,

.. code-block:: yaml

  ...
    data:
      -
        data_type: spaceranger
        data_path: /path/to/project_1/dataset_1/spaceranger/output/
        args:
          load_embeddings: True # for spaceranger
          compute_embeddings: False # for intermediate h5ad


.. _vitessce_options:

Vitessce options
================

The ``vitessce_options`` map is used to write Vitessce config files.
One Vitessce config file is generated per `dataset`_.
Include relevant information from your data to be visualized.
All values are optional as they depend on them existing in your data.

Values that can be specified are as follows

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
metadata within the AnnData object. It is written directly to the Vitessce
config file. If they're incorrectly specified then an error will occur when
Vitessce tries to load the data.

The output config file can be manually edited without re-running the pipeline
to fix or adapt the visualization to your needs.
