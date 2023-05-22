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

   wget https://github.com/haniffalab/webatlas-pipeline/archive/refs/tags/v0.3.2.tar.gz

.. code-block:: shell
   :caption: Output
    
   Resolving github.com (github.com)... 140.82.121.3
   Connecting to github.com (github.com)|140.82.121.3|:443... connected.
   HTTP request sent, awaiting response... 302 Found
   Location: https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.3.2 [following]
   --2023-05-18 09:30:15--  https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.3.2
   Resolving codeload.github.com (codeload.github.com)... 140.82.121.9
   Connecting to codeload.github.com (codeload.github.com)|140.82.121.9|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: unspecified [application/x-gzip]
   Saving to: ‘v0.3.2.tar.gz’

   v0.3.2.tar.gz                                                      [  <=>                                                                                                                                               ]   2.70M  9.12MB/s    in 0.3s    

   2023-05-18 09:30:16 (9.12 MB/s) - ‘v0.3.2.tar.gz’ saved [2835534]

**#3. Extract the WebAtlas Pipeline**

Download and extract a  of the WebAtlas repo and change directory into the new repo: 

.. code-block:: shell
   :caption: Input

   tar -xzvf ./v0.3.2.tar.gz
   cd webatlas-pipeline-0.3.2

.. code-block:: shell
   :caption: Output
    
   webatlas-pipeline-0.3.2/
   webatlas-pipeline-0.3.2/.github/
   ...
   ...
   webatlas-pipeline-0.3.2/tests/input/simple_config.json
   webatlas-pipeline-0.3.2/tests/test_class.py

**#4. Check java is installed**

Make sure java 11 later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   java -version

.. code-block:: shell
   :caption: Output
   
   openjdk version "11.0.18" 2023-01-17
   OpenJDK Runtime Environment (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1)
   OpenJDK 64-Bit Server VM (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1, mixed mode, sharing)

**#5. Install Nextflow**

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

**#6. Check Docker is installed**

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

**#7. Download the sample data**

.. code-block:: shell
   :caption: Input

   mkdir -p ./input/CytAssist_FFPE_Human_Breast_Cancer
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_image.tif -O ./input/CytAssist_FFPE_Human_Breast_Cancer/image.tif
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_tissue_image.tif -O ./input/CytAssist_FFPE_Human_Breast_Cancer/tissue_image.tif
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_possorted_genome_bam.bam.bai -O ./input/CytAssist_FFPE_Human_Breast_Cancer/possorted_genome_bam.bam.bai
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_analysis.tar.gz -O ./input/CytAssist_FFPE_Human_Breast_Cancer/analysis.tar.gz
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5 -O ./input/CytAssist_FFPE_Human_Breast_Cancer/filtered_feature_bc_matrix.h5
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_raw_feature_bc_matrix.h5 -O ./input/CytAssist_FFPE_Human_Breast_Cancer/raw_feature_bc_matrix.h5
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz -O ./input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.tar.gz -O ./input/CytAssist_FFPE_Human_Breast_Cancer/filtered_feature_bc_matrix.tar.gz
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_raw_feature_bc_matrix.tar.gz -O ./input/CytAssist_FFPE_Human_Breast_Cancer/raw_feature_bc_matrix.tar.gz
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_molecule_info.h5 -O ./input/CytAssist_FFPE_Human_Breast_Cancer/molecule_info.h5
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_cloupe.cloupe -O ./input/CytAssist_FFPE_Human_Breast_Cancer/cloupe.cloupe
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_possorted_genome_bam.bam -O ./input/CytAssist_FFPE_Human_Breast_Cancer/possorted_genome_bam.bam
   wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_metrics_summary.csv -O ./input/CytAssist_FFPE_Human_Breast_Cancer/metrics_summary.csv

.. code-block:: shell
   :caption: Output

   --2023-05-17 21:37:57--  https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_metrics_summary.csv
   Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 104.18.0.173, 104.18.1.173, 2606:4700::6812:ad, ...
   Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|104.18.0.173|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: 803 [text/csv]
   Saving to: ‘./input/CytAssist_FFPE_Human_Breast_Cancer/metrics_summary.csv’

   ./input/CytAssist_FFPE_Human_Breas 100%[================================================================>]     803  --.-KB/s    in 0s      

   2023-05-17 21:37:58 (7.16 MB/s) - ‘./input/CytAssist_FFPE_Human_Breast_Cancer/metrics_summary.csv’ saved [803/803]

**#8. Extract and process sample data**

.. code-block:: shell
   :caption: Input

   tar -xzvf ./input/CytAssist_FFPE_Human_Breast_Cancer/analysis.tar.gz -C ./input/CytAssist_FFPE_Human_Breast_Cancer
   tar -xzvf ./input/CytAssist_FFPE_Human_Breast_Cancer/spatial.tar.gz -C ./input/CytAssist_FFPE_Human_Breast_Cancer
   tar -xzvf ./input/CytAssist_FFPE_Human_Breast_Cancer/filtered_feature_bc_matrix.tar.gz -C ./input/CytAssist_FFPE_Human_Breast_Cancer
   tar -xzvf ./input/CytAssist_FFPE_Human_Breast_Cancer/raw_feature_bc_matrix.tar.gz -C ./input/CytAssist_FFPE_Human_Breast_Cancer
   cp ./input/CytAssist_FFPE_Human_Breast_Cancer/spatial/tissue_positions.csv ./input/CytAssist_FFPE_Human_Breast_Cancer/spatial/tissue_positions_list.csv

.. code-block:: shell
   :caption: Output

   analysis/umap/
   analysis/umap/gene_expression_2_components/
   ...
   ...
   raw_feature_bc_matrix/barcodes.tsv.gz
   raw_feature_bc_matrix/matrix.mtx.gz

**#9. Run the pipeline**

.. code-block:: shell
   :caption: Input

   NXF_VER=22.04.5 ./nextflow run main.nf -params-file templates/examples/CytAssist_FFPE_Human_Breast_Cancer.yaml -entry Full_pipeline

.. code-block:: shell
   :caption: Output

   N E X T F L O W  ~  version 22.04.5
   Launching `main.nf` [insane_dijkstra] DSL2 - revision: 1b6a73f4d6
   [05/d2276b] process > Full_pipeline:Process_files:route_file (spaceranger, CytAssist_FFPE_Human_Breast_Cancer)   [100%] 1 of 1 ✔
   [0c/3ffdac] process > Full_pipeline:Process_images:Generate_image ([visium, breast-cancer], label, CytAssist_... [100%] 1 of 1 ✔
   [f1/efaaae] process > Full_pipeline:Process_images:image_to_zarr (tissue_image.tif)                              [100%] 2 of 2 ✔
   [44/2bcaeb] process > Full_pipeline:Process_images:ome_zarr_metadata (METADATA.ome.xml)                          [100%] 2 of 2 ✔
   [43/04893d] process > Full_pipeline:Output_to_config:Build_config ([visium, breast-cancer])                      [100%] 1 of 1 ✔

   {"dimOrder": "XYZCT", "channel_names": [], "X": "19505", "Y": "21571", "Z": "1", "C": "1", "T": "1"}

   {"dimOrder": "XYZCT", "channel_names": [], "X": "19505", "Y": "21571", "Z": "1", "C": "3", "T": "1"}

**#10. Check execution was successful**

The output from the pipeline will indicate if the execution was successful. You can also
verify the expected directories are created. 

.. code-block:: shell
   :caption: Input

   ls -l ./output/CytAssist_FFPE_Human_Breast_Cancer/0.3.2

.. code-block:: shell
   :caption: Output

   total 1103476
   -rw-r--r--  1 ndh74 ndh74 288446018 May 17 21:42 tmp-visium-breast-cancer.h5ad
   drwxrwxr-x 11 ndh74 ndh74      4096 May 17 21:42 visium-breast-cancer-anndata.zarr
   -rw-r--r--  1 ndh74 ndh74      4667 May 17 21:43 visium-breast-cancer-config.json
   -rw-r--r--  1 ndh74 ndh74 841484966 May 17 21:42 visium-breast-cancer-label.tif
   drwxrwxr-x  4 ndh74 ndh74      4096 May 17 21:43 visium-breast-cancer-label.zarr
   drwxrwxr-x  4 ndh74 ndh74      4096 May 17 21:43 visium-breast-cancer-raw.zarr

**#11. Serve the data output through a local web server**

To browse and explore the data, you need to serve the output data through a web server.
You can use your preferred web server, but you must ensure the data is served over port 3000, 
at http://localhost:3000, and that CORS is enabled via the Access-Control-Allow-Origin header.

.. code-block:: shell
   :caption: Input

   npx http-server ./output/CytAssist_FFPE_Human_Breast_Cancer/0.3.2 --port 3000 --cors

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

**#12. Explore data in your browser**

Start your web browser and open:

https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=http://127.0.0.1:3000/visium-breast-cancer-config.json