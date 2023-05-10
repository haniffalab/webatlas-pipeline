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
    
**#3. Check java is installed**

Make sure java 11 later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   java -version
   openjdk version "11"

**#4. Install Nextflow**

Enter this command in your terminal to install nextflow in the current directory:

.. code-block:: shell
   :caption: Input

   curl -s https://get.nextflow.io | bash

**#5. Check Docker is installed**

Make sure Docker Engine 20.10 later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   docker version

.. code-block:: shell
   :caption: Output

   Client: Docker Engine - Community

**#6. Build the docker images**

.. code-block:: shell
   :caption: Input

   cd docker
   ./build-docker-imgs.sh
   cd ../

.. code-block:: shell
   :caption: Output

   somet

**#7. Download the data**

.. code-block:: shell
   :caption: Input

   wget https://www.something.com
   tar ...

.. code-block:: shell
   :caption: Output

   dgdf

**#8. Download the input configuration file**

.. code-block:: shell
   :caption: Input

   wget https://www.something.com
   tar ...

.. code-block:: shell
   :caption: Output

   dgdf

**#9. Run the pipeline**

.. code-block:: shell
   :caption: Input

   nextflow run main.nf -params-file visium.yaml -entry Full_pipeline

.. code-block:: shell
   :caption: Output

   dgdf

**#10. Check execution was successful**

**#11. Serve the data output through a web browser**

**#12. Load data in your browser**

