.. _multimodal_visualise:

Multimodal visualision
======================

The pipeline generates a Vitessce view config file for each ``project`` to visualise its multiple ``datasets`` together.
This file can then be used to load the views and data as configured in the parameters files.
Unlike the main conversion pipeline, the layout of the view config file generated for multimodal datasets cannot be configured.

For each spatial ``dataset`` (has ``is_spatial: true`` and contains either a ``raw_image`` or ``label_image``) the pipeline adds an spatial image component and a layer controller to the view config.
For each non-spatial ``dataset`` (has ``is_spatial: false`` or contains no images) the pipeline adds a scatterplot component to the view config to visualise an embedding.
As datasets have been subsetted to the shared features only, the pipeline adds one feature selection component to the view config from which the user can query across datasets.
If an additional feature has been specified, the pipeline adds another feature selection component to the view config for that particular feature type.

To set the layout of the view config the pipeline just concatenates components in a recursive manner.
Thus, components might end up with sizes that are not the most ideal.
During visualisation the user can resize and reorganise the components as needed.
Alternatively, the generated view config file can be manually modified to set a more appropriate layout.

You can locally serve and visualize the data in a few steps.

By default, the base ``url`` used within the view config is ``http://localhost:3000/`` 
(this can be changed in the :ref:`parameters file <multimodal_parameters_file>`).
This ``url`` indicates Vitessce to look for data at that location.

You can set up a ``http`` server locally to serve the processed files so a Vitessce instance can load them.
There are several tools that can setup a ``http`` server.
We recommend using `serve <https://www.npmjs.com/package/serve>`__ (requires `Node.js <https://nodejs.org/en/>`__),
but you can use any tool that can enable CORS. 

You can serve the view config file and data by specifying the output directory
(note that the pipeline adds its version to the ``outdir`` defined in the :ref:`parameters file <multimodal_parameters_file>` file). 

.. parsed-literal::

   serve -C -p 3000 /path/to/outdir/|release|/

Make sure to enable CORS and set the appropriate port number.
In this case, using `serve <https://www.npmjs.com/package/serve>`__, this is done through the ``-C`` and ``-p`` flags respectively.

Your view configs should then be accessible at ``http://localhost:3000/{project}-{dataset}-config.json``.

You can then load them in a Vitessce instance like the `WebAtlas app <https://github.com/haniffalab/webatlas-app>`__ 
deployed at `<https://webatlas.cog.sanger.ac.uk/latest/index.html>`__.

Specify your locally served view config through the ``config`` parameter like
``https://webatlas.cog.sanger.ac.uk/latest/index.html?config=http://localhost:3000/{project}-{dataset}-config.json``
and load this URL in your browser to visualize your data in a Vitessce viewer.