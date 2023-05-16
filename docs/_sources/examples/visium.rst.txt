.. _example_visium:

Visium
======

Follow this step by step guide to process and visualise a Visium sample in your web browser. 
It can be used on any POSIX compatible system (Linux, OS X, etc).

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

**#6. Build the docker images**

.. code-block:: shell
   :caption: Input

   cd docker
   ./build-docker-imgs.sh
   cd ../

.. code-block:: shell
   :caption: Output

   [+] Building 152.0s (9/9) FINISHED                                                                                                        
   => [internal] load build definition from Dockerfile                                                                                 0.0s
   => => transferring dockerfile: 874B                                                                                                 0.0s
   => [internal] load .dockerignore                                                                                                    0.0s
   => => transferring context: 2B                                                                                                      0.0s
   => [internal] load metadata for docker.io/library/python:3.9                                                                        1.7s
   => [1/4] FROM docker.io/library/python:3.9@sha256:6ea9dafc96d7914c5c1d199f1f0195c4e05cf017b10666ca84cb7ce8e2699d51                 26.9s
   => => resolve docker.io/library/python:3.9@sha256:6ea9dafc96d7914c5c1d199f1f0195c4e05cf017b10666ca84cb7ce8e2699d51                  0.0s
   => => sha256:4eedd9c5abf7e5f63753a5e788cb0872a715fa1141e8ce5ea87638e6cd370a41 54.58MB / 54.58MB                                    11.3s
   => => sha256:6ea9dafc96d7914c5c1d199f1f0195c4e05cf017b10666ca84cb7ce8e2699d51 1.86kB / 1.86kB                                       0.0s
   => => sha256:c65dadac8789fed40962578392e99a0528dcb868442c75d144e68ba858984837 2.01kB / 2.01kB                                       0.0s
   => => sha256:27ac39eccd1fd7d6c7f78eefc98233e633a7178dd05910983ba905ecf71561f3 7.51kB / 7.51kB                                       0.0s
   => => sha256:5d79063a01c561833dc6546d4e647fda0121a59e1a9a17874a3e30854555475e 15.76MB / 15.76MB                                     3.0s
   => => sha256:918547b9432687b1e1d238e82dc1e0ea0b736aafbf3c402eea98c6db81a9cb65 55.05MB / 55.05MB                                     6.5s
   => => sha256:9cdadd40055fb82fef74a0a38f29fb1ee7b6922e18370ebdc3502ec66f79fa3f 196.85MB / 196.85MB                                  22.0s
   => => sha256:2a12d0031f3fc457f41eaa0e29c0eaf4ba1376735aa940268ebccc76abb49f04 6.29MB / 6.29MB                                       8.7s
   => => extracting sha256:918547b9432687b1e1d238e82dc1e0ea0b736aafbf3c402eea98c6db81a9cb65                                            0.7s
   => => extracting sha256:5d79063a01c561833dc6546d4e647fda0121a59e1a9a17874a3e30854555475e                                            0.2s
   => => sha256:cea461a97d8717d5ef1eb563adaf205de17b9bf5b9bd57bc6bb2eae2296a8740 16.15MB / 16.15MB                                    11.9s
   => => extracting sha256:4eedd9c5abf7e5f63753a5e788cb0872a715fa1141e8ce5ea87638e6cd370a41                                            1.5s
   => => sha256:a48c72dfa8c4e0da323ce32dac703c79d8f2de134e6a420cdfbc868323cfc68a 244B / 244B                                          12.0s
   => => sha256:c343b921680a19bf22934b83c02f42d8fdd0d8559417f63c1dcd2d7997c22da3 2.90MB / 2.90MB                                      12.6s
   => => extracting sha256:9cdadd40055fb82fef74a0a38f29fb1ee7b6922e18370ebdc3502ec66f79fa3f                                            2.9s
   => => extracting sha256:2a12d0031f3fc457f41eaa0e29c0eaf4ba1376735aa940268ebccc76abb49f04                                            0.2s
   => => extracting sha256:cea461a97d8717d5ef1eb563adaf205de17b9bf5b9bd57bc6bb2eae2296a8740                                            0.3s
   => => extracting sha256:a48c72dfa8c4e0da323ce32dac703c79d8f2de134e6a420cdfbc868323cfc68a                                            0.0s
   => => extracting sha256:c343b921680a19bf22934b83c02f42d8fdd0d8559417f63c1dcd2d7997c22da3                                            0.1s
   => [internal] load build context                                                                                                    0.0s
   => => transferring context: 736B                                                                                                    0.0s
   => [2/4] COPY ./requirements.txt /requirements.txt                                                                                  0.8s
   => [3/4] RUN apt-get update &&     apt-get install -y --no-install-recommends       wget unzip cmake g++ ant procps libblosc1 lib  27.6s
   => [4/4] RUN pip install --upgrade pip setuptools distlib --no-cache-dir &&     pip install --no-cache-dir -r /requirements.txt    91.3s
   => exporting to image                                                                                                               3.7s
   => => exporting layers                                                                                                              3.7s
   => => writing image sha256:71bf0c89267114ae8df748dfd8d262c253f7cd2658c2221d329e50604cab1722                                         0.0s
   => => naming to docker.io/haniffalab/webatlas-pipeline:0.0.1 

**#7. Download and extract the sample data**

.. code-block:: shell
   :caption: Input

   wget https://www.something.com/WSSS_THYst9699525.tar.gz
   mkdir -p ./input/WSSS_THYst9699525
   tar -xzvf WSSS_THYst9699525.tar.gz -C ./input/WSSS_THYst9699525

.. code-block:: shell
   :caption: Output

   ./
   ./WSSS_THYst9699525.h5ad
   ./WSSS_THYst9699525.yaml
   ./S3364_C59-FLEG-FO1_C1_HnE-2020-10-07-10.08.32.tif

**#8. Run the pipeline**

.. code-block:: shell
   :caption: Input

   nextflow run main.nf -params-file input/WSSS_THYst9699525/WSSS_THYst9699525.yaml -entry Full_pipeline

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


**#9. Check execution was successful**

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

**#10. Serve the data output through a local web server**

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

**#11. Explore data in your browser**

Start your web browser and open:

https://webatlas.cog.sanger.ac.uk/dev/index.html?theme=dark&config=http://127.0.0.1:8080/lowerlimb-visium-WSSS_THYst9699525-config.json