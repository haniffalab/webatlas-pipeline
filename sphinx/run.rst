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

To run the ``Full_pipeline`` use

.. code-block:: sh

   nextflow run main.nf -params-file /path/to/params.yaml -entry Full_pipeline


This will handle all input files, whether they are data files or images, for all datasets
defined.

You can modify the entry point if you're interested in only getting the converted outputs.
Use ``-entry Process_files`` or ``-entry Process_images`` as you need.

Further reading
---------------

Docker image pulling/local conda env creation are handled by Nextflow.
Please refer to
`this <https://www.nextflow.io/docs/latest/getstarted.html>`__ for
detailed information.