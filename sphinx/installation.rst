.. _official nextflow documentation: https://www.nextflow.io/index.html#GetStarted
.. _official Docker Install guide: https://docs.docker.com/engine/install/
.. _releases on GitHub: https://github.com/haniffalab/webatlas-pipeline/releases
.. _conda: https://docs.conda.io/projects/miniconda/en/latest/
.. _mamba: https://mamba.readthedocs.io/en/latest/mamba-installation.html

.. _installation:

Installation
============

Download the WebAtlas Pipeline release. You can look for previous `releases on GitHub`_.

Using `wget`
""""""""""""

.. code-block:: shell
   :caption: Input

   wget https://github.com/haniffalab/webatlas-pipeline/archive/refs/tags/v0.5.3.tar.gz

.. code-block:: shell
   :caption: Expected Output

   Resolving github.com (github.com)... 140.82.121.3
   Connecting to github.com (github.com)|140.82.121.3|:443... connected.
   HTTP request sent, awaiting response... 302 Found
   Location: https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.5.3 [following]
   --2023-05-18 09:30:15--  https://codeload.github.com/haniffalab/webatlas-pipeline/tar.gz/refs/tags/v0.5.3
   Resolving codeload.github.com (codeload.github.com)... 140.82.121.9
   Connecting to codeload.github.com (codeload.github.com)|140.82.121.9|:443... connected.
   HTTP request sent, awaiting response... 200 OK
   Length: unspecified [application/x-gzip]
   Saving to: ‘v0.5.3.tar.gz’

   v0.5.3.tar.gz [ <=>                                           ]   2.70M  9.12MB/s    in 0.3s    

   2023-05-18 09:30:16 (9.12 MB/s) - ‘v0.5.3.tar.gz’ saved [2835534]


Using `curl`
""""""""""""

.. code-block:: shell
   :caption: Input

   curl -L -o v0.5.3.tar.gz https://github.com/haniffalab/webatlas-pipeline/archive/refs/tags/v0.5.3.tar.gz

.. code-block:: shell
   :caption: Expected Output

   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                    Dload  Upload   Total   Spent    Left  Speed
   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
   100 2844k    0 2844k    0     0  1970k      0 --:--:--  0:00:01 --:--:-- 2539k


Extract the WebAtlas compressed tag and change directory into the new repo.

.. code-block:: shell
   :caption: Input

   tar -xzvf ./v0.5.3.tar.gz
   cd webatlas-pipeline-0.5.3

.. code-block:: shell
   :caption: Expected Output
    
   webatlas-pipeline-0.5.3/
   webatlas-pipeline-0.5.3/.github/
   ...
   ...
   webatlas-pipeline-0.5.3/tests/input/simple_config.json
   webatlas-pipeline-0.5.3/tests/test_class.py

.. _environment:

Environment setup
=================

.. _environment_conda:

Follow these Environment Setup instructions using conda or manually installing the required components to run WebAtlas.

Using conda
-----------

If you have `conda`_ or `mamba`_ already installed then you can use the ``environment.yaml`` file included in the WebAtlas release to create the environment.

.. code-block:: shell
   :caption: Input

   conda env create -f envs/environment.yaml

Then make sure you activate the ``webatlas`` environment before you use the pipeline.

.. code-block:: shell
   :caption: Input

   conda activate webatlas

.. warning::
   Users working on newer Silicon-based Macs may encounter problems installing this environment.
   Some packages have not yet been compiled for Apple silicon processors therefore, 
   we recommend you install the packages originally compiled for Mac computers with Intel processors. Set
   an environment variable that specifies the architecture before installing and activating the Conda
   environment, like this:

   .. code-block:: shell
      :caption: Input

      export CONDA_SUBDIR=osx-64
      conda env create -f envs/environment.yaml 
      conda activate webatlas

.. _environment_manual:

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

**#2. Install Nextflow**

Java is required by Nextflow. Refer to `Nextflow's guidelines <https://www.nextflow.io/docs/latest/getstarted.html#requirements>`__ to install it.

Enter the following command in your terminal to install Nextflow in the current directory:

Using `wget`
""""""""""""

.. code-block:: shell
   :caption: Input

   wget -qO- https://get.nextflow.io | bash


Using `curl`
""""""""""""

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

**#3. Check Docker is installed (optional)**

If you want to use Docker, make sure Docker Engine 20.10 or later is installed on your computer by using the command:

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
