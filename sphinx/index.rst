|Tests| |Sphinx| |Coverage| |DOI|
 
.. |Tests| image:: https://github.com/haniffalab/webatlas-pipeline/actions/workflows/tests-python.yml/badge.svg
   :target: https://github.com/haniffalab/webatlas-pipeline/actions/workflows/tests-python.yml
.. |Sphinx| image:: https://github.com/haniffalab/webatlas-pipeline/actions/workflows/sphinx-build.yml/badge.svg
   :target: https://github.com/haniffalab/webatlas-pipeline/actions/workflows/sphinx-build.yml
.. |Coverage| image:: https://codecov.io/gh/haniffalab/webatlas-pipeline/branch/main/graph/badge.svg?token=7HQVFH08WJ
   :target: https://app.codecov.io/gh/haniffalab/webatlas-pipeline
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7405818.svg
   :target: https://doi.org/10.5281/zenodo.7405818

WebAtlas pipeline
=================

This Nextflow pipeline processes spatial and single-cell experiment data for visualisation in `WebAtlas App`_. 
The pipeline generates data files for `supported data types`_, and builds a `view config`_.

.. _WebAtlas App: https://github.com/haniffalab/webatlas-app
.. _supported data types: https://vitessce.io/docs/data-types-file-types/
.. _view config: https://vitessce.io/docs/view-config-json/

.. toctree::
   :maxdepth: 1
   :caption: Documentation
   :glob:

   installation
   configuration
   multimodal_configuration
   run
   visualise
   Demos <https://cellatlas.io/webatlas>
   testing
   modules

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Multimodal
  
   multimodal/overview
   multimodal/configuration
   multimodal/run
   multimodal/visualise

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Example Workflows
  
   examples/visium
   examples/xenium 

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Project Links
  
   citing
   Source Code <https://github.com/haniffalab/webatlas-pipeline>
   Issue Tracker <https://github.com/haniffalab/webatlas-pipeline/issues>
   WebAtlas App <https://github.com/haniffalab/webatlas-app>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Project Teams
  
   Bayraktar Lab <https://www.sanger.ac.uk/group/bayraktar-group>
   Haniffa Lab <https://haniffalab.com>
   Open Microscopy Environment <https://www.openmicroscopy.org/>
