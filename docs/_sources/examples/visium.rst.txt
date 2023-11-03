.. _example_visium:

Visium
======

Sample details
**************

.. list-table::
   :widths: 25 75
   :header-rows: 0

   * - Study Name
     - `High resolution mapping of the breast cancer tumor microenvironment using integrated single cell, spatial and in situ analysis of FFPE tissue <https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast>`__
   * - WebAtlas
     - `Demo <https://webatlas.cog.sanger.ac.uk/latest/index.html?config=https://webatlas.cog.sanger.ac.uk/configs/dev/visium/human/breast/cancer/config.json>`__     
   * - Tissue
     - Human breast cancer
   * - Data Source Link
     - `CytAssist_FFPE_Human_Breast_Cancer <https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast>`__

Steps to reproduce
******************

Follow the steps below to reproduce this sample in the pipeline, and visualise the data yourself 
in your web browser. It can be followed on any POSIX compatible system (Linux, OS X, etc). This
example requires you to have already :ref:`setup your environment first <environment>`.

**#1. Download the sample data**

.. code-block:: shell
   :caption: Input

   mkdir -p input/CytAssist_FFPE_Human_Breast_Cancer
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_tissue_image.tif -O input/CytAssist_FFPE_Human_Breast_Cancer/tissue_image.tif
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz -O input/CytAssist_FFPE_Human_Breast_Cancer/analysis.tar.gz
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5 -O input/CytAssist_FFPE_Human_Breast_Cancer/filtered_feature_bc_matrix.h5
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz -O input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz

.. code-block:: shell
   :caption: Output

   --2023-05-17 21:37:57--  https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz -O ./input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz
   Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 104.18.0.173, 104.18.1.173, 2606:4700::6812:ad, ...
   Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|104.18.0.173|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: 34479952 (33M) [application/x-tar]
   Saving to: ‘input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz

   ./input/CytAssist_FFPE_Human_Breas 100%[================================================================>]  32.88M  --.-KB/s    in 0s      

   2023-05-17 21:37:58 (7.16 MB/s) - ‘input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz’ saved [34479952/34479952]

**#2. Extract and process sample data**

.. code-block:: shell
   :caption: Input

   tar -xzvf input/CytAssist_FFPE_Human_Breast_Cancer/analysis.tar.gz -C input/CytAssist_FFPE_Human_Breast_Cancer
   tar -xzvf input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz -C input/CytAssist_FFPE_Human_Breast_Cancer

.. code-block:: shell
   :caption: Output

   analysis/
   analysis/umap/
   analysis/umap/gene_expression_2_components/
   analysis/umap/gene_expression_2_components/projection.csv
   ...
   spatial/scalefactors_json.json
   spatial/aligned_fiducials.jpg
   spatial/tissue_hires_image.png

**#3. Run the pipeline**

.. code-block:: shell
   :caption: Input

   nextflow run main.nf \
         -params-file templates/examples/CytAssist_FFPE_Human_Breast_Cancer.yaml \
         -entry Full_pipeline

.. code-block:: shell
   :caption: Output

   N E X T F L O W  ~  version 22.04.5
   Launching `main.nf` [insane_dijkstra] DSL2 - revision: 1b6a73f4d6
   [05/d2276b] process > Full_pipeline:Process_files:route_file (spaceranger, CytAssist_FFPE_Human_Breast_Cancer)   [100%] 1 of 1 ✔
   [0c/3ffdac] process > Full_pipeline:Process_images:Generate_image ([visium, breast-cancer], label, CytAssist_... [100%] 1 of 1 ✔
   [f1/efaaae] process > Full_pipeline:Process_images:image_to_zarr (tissue_image.tif)                              [100%] 2 of 2 ✔
   [44/2bcaeb] process > Full_pipeline:Process_images:ome_zarr_metadata (METADATA.ome.xml)                          [100%] 2 of 2 ✔
   [43/04893d] process > Full_pipeline:Output_to_config:Build_config ([visium, breast-cancer])                      [100%] 1 of 1 ✔

**#4. Check execution was successful**

The output from the pipeline will indicate if the execution was successful. You can also
verify the expected directories are created. 

.. code-block:: shell
   :caption: Input

   ls -l output/CytAssist_FFPE_Human_Breast_Cancer/0.4.0

.. code-block:: shell
   :caption: Output

   total 1103476
   -rw-r--r--  1 dh74 dh74 288446018 May 17 21:42 tmp-visium-breast-cancer.h5ad
   drwxrwxr-x 11 dh74 dh74      4096 May 17 21:42 visium-breast-cancer-anndata.zarr
   -rw-r--r--  1 dh74 dh74      4667 May 17 21:43 visium-breast-cancer-config.json
   -rw-r--r--  1 dh74 dh74 841484966 May 17 21:42 visium-breast-cancer-label.tif
   drwxrwxr-x  4 dh74 dh74      4096 May 17 21:43 visium-breast-cancer-label.zarr
   drwxrwxr-x  4 dh74 dh74      4096 May 17 21:43 visium-breast-cancer-raw.zarr

**#5. Serve the data output through a local web server**

To browse and explore the data, you need to serve the output data through a web server.
You can use your preferred web server, but you must ensure the data is served over port 3000, 
at http://localhost:3000, and that CORS is enabled via the Access-Control-Allow-Origin header.

.. code-block:: shell
   :caption: Input

   npx http-server output/CytAssist_FFPE_Human_Breast_Cancer/0.4.0 --port 3000 --cors

.. code-block:: shell
   :caption: Output

   Starting up http-server, serving ./

   http-server version: 14.1.1

   http-server settings: 
   CORS: true
   Cache: 3600 seconds
   Connection Timeout: 120 seconds
   Directory Listings: visible
   AutoIndex: visible
   Serve GZIP Files: false
   Serve Brotli Files: false
   Default File Extension: none

   Available on:
   http://127.0.0.1:3000
   http://192.168.0.23:3000
   Hit CTRL-C to stop the server

**#6. Explore data in your browser**

Start your web browser and open:

https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=http://127.0.0.1:3000/visium-breast-cancer-config.json