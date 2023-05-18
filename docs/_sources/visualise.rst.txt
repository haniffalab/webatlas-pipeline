.. _visualise:

Visualising
===========

The pipeline generates a Vitessce view config file for each processed dataset.
This file can then be used to load the views and data as configured in the parameters files.

You can locally serve and visualize the data in a few steps.

By default, the base ``url`` used within the view config is ``http://localhost:3000/`` 
(this can be changed in the :ref:`parameters file <parameters_file>`).
This ``url`` indicates Vitessce to look for data at that location.

You can set up an http server locally to serve the processed files so a Vitessce instance can load them.
There are several tools that can setup an http server.
We recommend using `serve <https://www.npmjs.com/package/serve>`__ (requires `Node.js <https://nodejs.org/en/>`__),
but you can use any tool that can enable CORS.

You can serve the view config file and data by specifying the output directory
(note that the pipeline adds its version to the ``outdir`` defined in the :ref:`parameters file <parameters_file>` file). 

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