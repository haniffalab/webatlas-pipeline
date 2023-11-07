.. _run:

Running
=======

The pipeline offers a workflow that processes files, images, and 
builds a Vitessce config file from the generated files.
Alternatively, the workflow to process files only, or the workflow to process images only  
can be called independently.

Each of these workflows work as entry points that can be specified when running the
pipeline through the command line.

- The ``Full_pipeline`` workflow runs the other workflows as needed and
  builds a Vitessce config file per dataset.
- The ``Process_files`` workflow handles data files and their conversions.
- The ``Process_images`` workflow handles image files and/or label image data and their conversions.

Configurations and data are input through a :ref:`parameters yaml file <parameters_file>`.
If specified in the parameters file, the pipeline can additionally write datasets to `SpatialData <https://spatialdata.scverse.org/en/latest/index.html>`_.

To run the ``Full_pipeline`` use

.. code-block:: shell

   nextflow run main.nf -params-file /path/to/params.yaml -entry Full_pipeline


This will handle all input files, whether they are data files or images, for all datasets
defined.

You can modify the entry point if you're interested in only getting the converted outputs.
Use ``-entry Process_files`` or ``-entry Process_images`` as you need.

Multimodal integration
----------------------

Additional to the main conversion pipeline, we offer a subsequent pipeline to process multiple datasets with matching features to be able to visualise and query across all of them in a single portal.
This pipeline can also process extra features (e.g. celltypes) to visualise them across datasets in addition to their expression matrices.

Configurations and data are input through a :ref:`parameters yaml file <multimodal_configuration>` (slightly different from the parameters file required by the main pipeline.)

To run this pipeline use

.. code-block:: shell

   nextflow run multimodal.nf -params-file /path/to/multimodal-params.yaml

Running using Docker 
--------------------

The default pipeline will run on local executor without any type of environment creation. To run the pipeline using Docker containers use the ``-profile docker`` option:

.. code-block:: shell

   nextflow run main.nf \
            -params-file /path/to/params.yaml \
            -entry Full_pipeline \
            -profile docker

Pulling the containers when the pipline is launched may take a few minutes.

Running using Singularity 
-------------------------

The default pipeline will run on local executor without any type of environment creation. To run the pipeline using Singularity containers use the ``-profile singularity`` option:

.. code-block:: shell

   nextflow run main.nf \
            -params-file /path/to/params.yaml \
            -entry Full_pipeline \
            -profile singularity

Pulling the containers when the pipline is launched may take a few minutes.

Running using Conda 
-------------------

The default pipeline will run on local executor without any type of environment creation. If you've already setup your conda environment you don't have to do anything else.

However, if you are working on a compute cluster you will need to make sure the conda environment is avaiable and active in your worker nodes. To run the pipeline using a new conda environment use the ``-profile conda`` option:

.. code-block:: shell

   nextflow run main.nf \
            -params-file /path/to/params.yaml \
            -entry Full_pipeline \
            -profile conda

Creating the environment when the pipleine is launched may take a few minutes.

Further reading
---------------

For more information about Docker image pulling/local conda env creation in Nextflow please refer to Nextflow's official docs for `containers <https://www.nextflow.io/docs/latest/container.html>`__ and `conda <https://www.nextflow.io/docs/latest/conda.html>`__.