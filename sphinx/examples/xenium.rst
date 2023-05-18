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

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto; margin-bottom: 20px;">
        <iframe src="https://www.youtube.com/embed/ScMzIvxBSi4" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>

Follow the steps below to reproduce this sample in the pipeline, and visualise the data yourself in your web browser. 
It can be followed on any POSIX compatible system (Linux, OS X, etc).

**#1. Check git is installed**

Make sure git 2.17 or later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   git --version

.. code-block:: shell
   :caption: Output

   git version 2.25.1

**#2. Download lthe WebAtlas Pipeline**

Download the WebAtlas Pipeline: 

.. code-block:: shell
   :caption: Input

   wget https://github.com/haniffalab/webatlas-pipeline/archive/refs/tags/v0.3.0.tar.gz

.. code-block:: shell
   :caption: Output
    
   Resolving github.com (github.com)... 140.82.121.3
   Connecting to github.com (github.com)|140.82.121.3|:443... connected.
   HTTP request sent, awaiting response... 302 Found
   Location: https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.3.0 [following]
   --2023-05-18 09:30:15--  https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.3.0
   Resolving codeload.github.com (codeload.github.com)... 140.82.121.9
   Connecting to codeload.github.com (codeload.github.com)|140.82.121.9|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: unspecified [application/x-gzip]
   Saving to: ‘v0.3.0.tar.gz’

   v0.3.0.tar.gz                                                      [  <=>                                                                                                                                               ]   2.70M  9.12MB/s    in 0.3s    

   2023-05-18 09:30:16 (9.12 MB/s) - ‘v0.3.0.tar.gz’ saved [2835534]

**#2. Extract the WebAtlas Pipeline**

Download and extract a  of the WebAtlas repo and change directory into the new repo: 

.. code-block:: shell
   :caption: Input

   tar -xzvf ./v0.3.0.tar.gz
   cd webatlas-pipeline-0.3.0

.. code-block:: shell
   :caption: Output
    
   webatlas-pipeline-0.3.0/
   webatlas-pipeline-0.3.0/.github/
   ...
   ...
   webatlas-pipeline-0.3.0/tests/input/simple_config.json
   webatlas-pipeline-0.3.0/tests/test_class.py

**#3. Check java is installed**

Make sure java 11 later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   java -version

.. code-block:: shell
   :caption: Output
   
   openjdk version "11.0.18" 2023-01-17
   OpenJDK Runtime Environment (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1)
   OpenJDK 64-Bit Server VM (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1, mixed mode, sharing)

**#4. Install Nextflow**

Enter this command in your terminal to install nextflow in the current directory:

.. code-block:: shell
   :caption: Input

   curl -s https://get.nextflow.io | bash

.. code-block:: shell
   :caption: Output
   
   CAPSULE: Downloading dependency org.apache.ivy:ivy:jar:2.5.1
   ...
   CAPSULE: Downloading dependency io.nextflow:nf-commons:jar:23.04.1
                                                                        
         N E X T F L O W
         version 23.04.1 build 5866
         created 15-04-2023 06:51 UTC (07:51 BST)
         cite doi:10.1038/nbt.3820
         http://nextflow.io


   Nextflow installation completed. Please note:
   - the executable file `nextflow` has been created in the folder: ./webatlas-pipeline
   - you may complete the installation by moving it to a directory in your $PATH

**#5. Check Docker is installed**

Make sure Docker Engine 20.10 later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   docker version

.. code-block:: shell
   :caption: Output

   Client: Docker Engine - Community
   Version:           23.0.4
   API version:       1.42
   Go version:        go1.19.8
   Git commit:        f480fb1
   Built:             Fri Apr 14 10:32:23 2023
   OS/Arch:           linux/amd64
   Context:           default

**#6. Download the sample data**

.. code-block:: shell
   :caption: Input

   mkdir -p ./input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs
   wget https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip -P ./input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs

.. code-block:: shell
   :caption: Output

   --2023-05-17 15:05:24--  https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip
   Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 104.18.0.173, 104.18.1.173, 2606:4700::6812:ad, ...
   Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|104.18.0.173|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: 9861155708 (9.2G) [application/zip]
   Saving to: ‘./input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip’

   Xenium_FFPE_Human_Breast_Cancer_Rep1 100%[===================================================================>]   9.18G  14.3MB/s    in 10m 6s  

   2023-05-17 15:15:31 (15.5 MB/s) - ‘./input/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip’ saved [9861155708/9861155708]

**#6. Extract the sample data**

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

**#7. Run the pipeline**

.. code-block:: shell
   :caption: Input

   NXF_VER=22.04.5 ./nextflow run main.nf -params-file templates/examples/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.yaml -entry Full_pipeline

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

**#8. Check execution was successful**

The output from the pipeline will indicate if the execution was successful. You can also
verify the expected directories are created. 

.. code-block:: shell
   :caption: Input

   ls -l ./output/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/0.3.0

.. code-block:: shell
   :caption: Output

    total 3566252
    drwxrwxr-x 11 ndh74 ndh74       4096 May 17 16:08 xenium-breast-cancer-anndata.zarr
    -rw-r--r--  1 ndh74 ndh74       4984 May 17 16:40 xenium-breast-cancer-config.json
    -rw-r--r--  1 ndh74 ndh74 3651814848 May 17 16:12 xenium-breast-cancer-label.tif
    drwxrwxr-x  4 ndh74 ndh74       4096 May 17 16:13 xenium-breast-cancer-label.zarr
    drwxrwxr-x  4 ndh74 ndh74       4096 May 17 16:40 xenium-breast-cancer-raw.zarr

**#9. Serve the data output through a local web server**

To browse and explore the data, you need to serve the output data through a web server.
You can use your preferred web server, but you must ensure the data is served over port 3000, 
at http://localhost:3000, and that CORS is enabled via the Access-Control-Allow-Origin header.

.. code-block:: shell
   :caption: Input

   npx http-server ./output/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/0.3.0 --port 3000 --cors

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