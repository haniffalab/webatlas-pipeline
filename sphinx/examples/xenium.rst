.. _example_xenium:

Xenium
======

Sample details
**************

.. list-table::
   :widths: 25 75
   :header-rows: 0

   * - Study Name
     - `High resolution mapping of the breast cancer tumor microenvironment using integrated single cell, spatial and in situ analysis of FFPE tissue <https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast>`__
   * - WebAtlas
     - `Demo <https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=https://webatlas.cog.sanger.ac.uk/configs/dev/xenium/human/breast/cancer/config.json>`__     
   * - Tissue
     - Human breast cancer
   * - Data Source Link
     - `Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip <https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast>`__

Steps to reproduce
******************

Follow the steps below to reproduce this sample in the pipeline, and visualise the data yourself 
in your web browser. It can be followed on any POSIX compatible system (Linux, OS X, etc). This
example requires you to have already :ref:`setup your environment first <environment>`.

**#1. Download the sample data**

.. code-block:: shell
   :caption: Input

   mkdir -p input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs
   wget https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip -P input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs

.. code-block:: shell
   :caption: Output

   --2023-05-17 15:05:24--  https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip
   Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 104.18.0.173, 104.18.1.173, 2606:4700::6812:ad, ...
   Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|104.18.0.173|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: 9861155708 (9.2G) [application/zip]
   Saving to: ‘input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip’

   Xenium_FFPE_Human_Breast_Cancer_Rep1 100%[===================================================================>]   9.18G  14.3MB/s    in 10m 6s  

   2023-05-17 15:15:31 (15.5 MB/s) - ‘input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip’ saved [9861155708/9861155708]

**#2. Extract the sample data**

.. code-block:: shell
   :caption: Input

   unzip input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip -d input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs

.. code-block:: shell
   :caption: Output

   Archive:  input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip
      creating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/
     inflating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/experiment.xenium  
      creating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/cell_feature_matrix/
     inflating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/cell_feature_matrix/barcodes.tsv.gz  
            ...
            ... 
     inflating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/metrics_summary.csv  
     inflating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/gene_panel.json  
     inflating: input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs/analysis_summary.html 

**#3. Run the pipeline**

.. code-block:: shell
   :caption: Input

   nextflow run main.nf \
         -params-file templates/examples/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.yaml \
         -entry Full_pipeline

.. code-block:: shell
   :caption: Output

    N E X T F L O W  ~  version 22.10.6
    Launching `main.nf` [gigantic_murdock] DSL2 - revision: 1b6a73f4d6
    [fc/782a3f] process > Full_pipeline:Process_files:route_file (xenium, outs)                              [100%] 1 of 1 ✔
    [b0/f5ff27] process > Full_pipeline:Process_images:Generate_image ([xenium, breast-cancer], label, outs) [100%] 1 of 1 ✔
    [2b/054048] process > Full_pipeline:Process_images:image_to_zarr (morphology.ome.tif)                    [100%] 2 of 2 ✔
    [07/5e37c4] process > Full_pipeline:Process_images:ome_zarr_metadata (METADATA.ome.xml)                  [100%] 2 of 2 ✔
    [c8/f2378c] process > Full_pipeline:Output_to_config:Build_config ([xenium, breast-cancer])              [100%] 1 of 1 ✔

    Completed at: 17-May-2023 16:40:58
    Duration    : 32m 47s
    CPU hours   : 0.6
    Succeeded   : 7

**#4. Check execution was successful**

The output from the pipeline will indicate if the execution was successful. You can also
verify the expected directories are created. 

.. code-block:: shell
   :caption: Input

   ls -l output/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/0.4.0

.. code-block:: shell
   :caption: Output

    total 3566252
    drwxrwxr-x 11 dh74 dh74       4096 May 17 16:08 xenium-breast-cancer-anndata.zarr
    -rw-r--r--  1 dh74 dh74       4984 May 17 16:40 xenium-breast-cancer-config.json
    -rw-r--r--  1 dh74 dh74 3651814848 May 17 16:12 xenium-breast-cancer-label.tif
    drwxrwxr-x  4 dh74 dh74       4096 May 17 16:13 xenium-breast-cancer-label.zarr
    drwxrwxr-x  4 dh74 dh74       4096 May 17 16:40 xenium-breast-cancer-raw.zarr

**#5. Serve the data output through a local web server**

To browse and explore the data, you need to serve the output data through a web server.
You can use your preferred web server, but you must ensure the data is served over port 3000, 
at http://localhost:3000, and that CORS is enabled via the Access-Control-Allow-Origin header.

.. code-block:: shell
   :caption: Input

   npx http-server output/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/0.4.0 --port 3000 --cors

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

**#10. Explore data in your browser**

Start your web browser and open:

https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=http://127.0.0.1:3000/xenium-breast-cancer-config.json