.. _example_visium:

Visium
======

Sample details
**************

.. list-table::
   :widths: 25 75
   :header-rows: 0

   * - Study Name
     - `A human embryonic limb cell atlas resolved in space and time <https://doi.org/10.1101/2022.04.27.489800>`__
   * - WebAtlas
     - `Demo Link <https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=https://webatlas.cog.sanger.ac.uk/configs/dev/lower-limb/modalities/visium/slide8/config.json>`__     
   * - Tissue
     - Developing Human Lower Limb   
   * - Data Source Link
     - `E-MTAB-10367 <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10367>`__

Steps to reproduce
******************

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto; margin-bottom: 20px;">
        <iframe src="https://www.youtube.com/embed/ScMzIvxBSi4" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>

Follow the steps below to reproduce this sample in the pipeline, and visualise the data yourself in your web browser. 
It can be followed on any POSIX compatible system (Linux, OS X, etc).

**#1. Check Git is installed**

Make sure git 2.17 or later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   git --version

.. code-block:: shell
   :caption: Output

   git version 2.25.1

**#2. Clone the WebAtlas pipeline repository**

Make a copy of the WebAtlas repo and change directory into the new repo: 

.. code-block:: shell
   :caption: Input

   git clone git@github.com:haniffalab/webatlas-pipeline.git
   cd webatlas-pipeline

.. code-block:: shell
   :caption: Output
    
   Cloning into 'webatlas-pipeline'...
   remote: Enumerating objects: 1966, done.
   remote: Counting objects: 100% (341/341), done.
   remote: Compressing objects: 100% (144/144), done.
   remote: Total 1966 (delta 227), reused 204 (delta 195), pack-reused 1625
   Receiving objects: 100% (1966/1966), 5.78 MiB | 5.91 MiB/s, done.
   Resolving deltas: 100% (1292/1292), done.

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

**#6. Download and extract the sample data**

.. code-block:: shell
   :caption: Input

   mkdir -p ./input/WSSS_THYst9699525
   wget https://storage.googleapis.com/haniffa-lab/resources/docs/webatlas/example-workflows/visium/WSSS_THYst9699525.tar.gz -P ./input/WSSS_THYst9699525

.. code-block:: shell
   :caption: Output

   --2023-05-17 12:43:58--  https://storage.googleapis.com/haniffa-lab/resources/docs/webatlas/example-workflows/visium/WSSS_THYst9699525.tar.gz
   Resolving storage.googleapis.com (storage.googleapis.com)... 172.217.169.16, 216.58.212.208, 216.58.212.240, ...
   Connecting to storage.googleapis.com (storage.googleapis.com)|172.217.169.16|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: 245291993 (234M) [application/gzip]
   Saving to: ‘WSSS_THYst9699525.tar.gz’

   WSSS_THYst9699525.tar.gz                  100%[=====================================================================================>] 233.93M  16.1MB/s    in 15s     

   2023-05-17 12:44:13 (15.8 MB/s) - ‘WSSS_THYst9699525.tar.gz’ saved [245291993/245291993]

.. code-block:: shell
   :caption: Input

   tar -xzvf ./input/WSSS_THYst9699525/WSSS_THYst9699525.tar.gz -C ./input/WSSS_THYst9699525

.. code-block:: shell
   :caption: Output

   ./
   ./WSSS_THYst9699525.h5ad
   ./WSSS_THYst9699525.yaml
   ./S3364_C59-FLEG-FO1_C1_HnE-2020-10-07-10.08.32.tif

**#7. Run the pipeline**

.. code-block:: shell
   :caption: Input

   ./nextflow run main.nf -params-file input/WSSS_THYst9699525/WSSS_THYst9699525.yaml -entry Full_pipeline

.. code-block:: shell
   :caption: Output

   N E X T F L O W  ~  version 22.10.6
   Launching `main.nf` [intergalactic_tuckerman] DSL2 - revision: 903c9797fa
   executor >  local (7)
   executor >  local (7)
   [42/8000d9] process > Full_pipeline:Process_files:route_file (h5ad, WSSS_THYst9699525.h5ad)                    [100%] 1 of 1 ✔
   [fa/e7a593] process > Full_pipeline:Process_images:Generate_image ([lowerlimb, visium-WSSS_THYst9699525], l... [100%] 1 of 1 ✔
   [42/2a6dab] process > Full_pipeline:Process_images:image_to_zarr (S3364_C59-FLEG-FO1_C1_HnE-2020-10-07-10.0... [100%] 2 of 2 ✔
   [4e/59100f] process > Full_pipeline:Process_images:ome_zarr_metadata (METADATA.ome.xml)                        [100%] 2 of 2 ✔
   [58/0f5162] process > Full_pipeline:Output_to_config:Build_config ([lowerlimb, visium-WSSS_THYst9699525])      [100%] 1 of 1 ✔

   lowerlimb-visium-WSSS_THYst9699525-anndata.zarr


   {"dimOrder": "XYZCT", "channel_names": [], "X": "15040", "Y": "26680", "Z": "1", "C": "1", "T": "1"}


   {"dimOrder": "XYZCT", "channel_names": [], "X": "15040", "Y": "26680", "Z": "1", "C": "3", "T": "1"}


**#8. Check execution was successful**

The output from the pipeline will indicate if the execution was successful. You can also
verify the expected directories are created. 

.. code-block:: shell
   :caption: Input

   ls -l ./output/WSSS_THYst9699525/0.0.1

.. code-block:: shell
   :caption: Output

   drwxrwxr-x 11 dave dave      4096 May 10 11:35 lowerlimb-WSSS_THYst9699525-anndata.zarr
   -rw-r--r--  1 dave dave      5104 May 10 11:36 lowerlimb-WSSS_THYst9699525-config.json
   -rw-r--r--  1 dave dave 802534656 May 10 11:35 lowerlimb-WSSS_THYst9699525-label.tif
   drwxrwxr-x  4 dave dave      4096 May 10 11:35 lowerlimb-WSSS_THYst9699525-label.zarr
   drwxrwxr-x  4 dave dave      4096 May 10 11:36 lowerlimb-WSSS_THYst9699525-raw.zarr

**#9. Serve the data output through a local web server**

To browse and explore the data, you need to serve the output data through a web server.
You can use your preferred web server, but you must ensure the data is served over port 3000, 
at http://localhost:3000, and that CORS is enabled via the Access-Control-Allow-Origin header.

.. code-block:: shell
   :caption: Input

   npx http-server ./output/WSSS_THYst9699525/0.0.1 --port 3000 --cors

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

https://webatlas.cog.sanger.ac.uk/latest/index.html?theme=dark&config=http://127.0.0.1:3000/lowerlimb-WSSS_THYst9699525-config.json