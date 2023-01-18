.. _run:

Running
=======

The pipeline contains workflows to process files (h5ad files,
spaceranger output, csv/tsv files), convert images to Zarrs and build a
Vitessce config file from the generated files.

-  The ``scRNAseq_pipeline`` entry point handles datasets with no image
   data. This processes files and builds a Vitessce config file.
-  The ``ISS_pipeline`` entry point handles datasets that contain both
   raw and label images.
-  The ``Visium_pipeline`` entry point handles datasets that contain a
   raw image and generates a label image from the specified h5ad file
   where it expects to have a ``spatial`` key within ``uns``.

Each dataset in a tsv file should belong to the same modality so they
can all be run through the corresponding pipeline.

::

   nextflow run main.nf -params-file [your_params].yaml -entry [modality]_pipeline

You can also just convert the images to Zarrs:

::

   nextflow run main.nf -params-file [your_params].yaml -entry To_ZARR

Or just process files:

::

   nextflow run main.nf -params-file [your_params].yaml -entry Process_files

Additionally, you can build a config file without processing files nor
images, using a slightly different `yaml
file <templates/config_template.yaml>`__, per dataset, where you would
provide filenames and necessary metadata:

::

   nextflow run main.nf -params-file [your_params].yaml -entry Config

Further reading:
----------------

Docker image pulling/local conda env creation are handled by nextflow.
Please refer to
`this <https://www.nextflow.io/docs/latest/getstarted.html>`__ for
detailed information.
