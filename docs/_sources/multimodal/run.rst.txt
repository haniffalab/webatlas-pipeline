.. _multimodal_run:

Multimodal run
===========

In additional to the main conversion pipeline, we offer a subsequent pipeline to process multiple datasets with matching features. This allows the users to 
visualise and query all common features such as genes and cell types across all modalities from a single web portal.

Configurations and data are input through a :ref:`parameters yaml file <multimodal_configuration>` (slightly different from the parameters file required by the main pipeline).

To run this pipeline use

.. code-block:: shell

   nextflow run multimodal.nf -params-file /path/to/multimodal-params.yaml

Running using Docker 
--------------------

The default pipeline will run on local executor without any type of environment creation. To run the pipeline using Docker containers use the ``-profile docker`` option:

.. code-block:: shell

   nextflow run multimodal.nf \
            -params-file /path/to/multimodal-params.yaml \
            -profile docker

Pulling the containers when the pipline is launched may take a few minutes.

Running using Singularity 
-------------------------

The default pipeline will run on local executor without any type of environment creation. To run the pipeline using Singularity containers use the ``-profile singularity`` option:

.. code-block:: shell

   nextflow run multimodal.nf \
            -params-file /path/to/multimodal-params.yaml \
            -profile singularity

Pulling the containers when the pipline is launched may take a few minutes.

Running using Conda 
-------------------

The default pipeline will run on local executor without any type of environment creation. If you've already setup your conda environment you don't have to do anything else.

However, if you are working on a compute cluster you will need to make sure the conda environment is avaiable and active in your worker nodes. To run the pipeline using a new conda environment use the ``-profile conda`` option:

.. code-block:: shell

   nextflow run multimodal.nf \
            -params-file /path/to/multimodal-params.yaml \
            -profile conda

Creating the environment when the pipleine is launched may take a few minutes.

Further reading
---------------

For more information about Docker image pulling/local conda env creation in Nextflow please refer to Nextflow's official docs for `containers <https://www.nextflow.io/docs/latest/container.html>`__ and `conda <https://www.nextflow.io/docs/latest/conda.html>`__.