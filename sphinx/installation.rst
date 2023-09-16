.. _official nextflow documentation: https://www.nextflow.io/index.html#GetStarted
.. _official Docker Install guide: https://docs.docker.com/engine/install/
.. _releases on GitHub: https://github.com/haniffalab/webatlas-pipeline/releases
.. _conda: https://docs.conda.io/projects/miniconda/en/latest/
.. _mamba: https://mamba.readthedocs.io/en/latest/mamba-installation.html

.. _installation:

Installation
============

Manual setup
------------

**#1. Check git is installed**

Make sure git 2.17 or later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   git --version

.. code-block:: shell
   :caption: Output

   git version 2.25.1

If Git is missing you will have to follow the `Getting Started Installing Git guide <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__.

**#2. Check java is installed**

Make sure Java 11 or later is installed on your computer by using the command:

.. code-block:: shell
   :caption: Input

   java -version

.. code-block:: shell
   :caption: Output

   openjdk version "11.0.18" 2023-01-17
   OpenJDK Runtime Environment (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1)
   OpenJDK 64-Bit Server VM (build 11.0.18+10-post-Ubuntu-0ubuntu120.04.1, mixed mode, sharing)

**#3. Install Nextflow**

Enter the following command in your terminal to install nextflow in the current directory:

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

You can read more about how to install nextflow in the `official nextflow documentation`_.


**#4. Check Docker is installed**

Make sure Docker Engine 20.10 or later is installed on your computer by using the command:

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

Follow the `official Docker Install guide`_ if it is not installed already.

**#5. Clone the repository**

This will clone the repository to your home directory under folder .nextflow/assets/:

.. code-block:: shell
   :caption: Input

   nextflow pull haniffalab/webatlas-pipeline

.. code-block:: shell
   :caption: Output

   Checking haniffalab/webatlas-pipeline ...
        downloaded from https://github.com/haniffalab/webatlas-pipeline.git - revision: 1131fcc69c [main]

**#6. Build local docker images (optional)**

When using docker the pipleine can use local images or pull them from DockerHub. If you want to build the images yourself you can do it like this:

::

    cd ~/.nextflow/assets/haniffalab/webatlas-pipeline/docker/
    ./build-docker-imgs.sh
