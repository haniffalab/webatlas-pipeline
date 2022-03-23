# sci-spatial-data

Utils to pre-process data for [Vitessce](http://github.com/hms-dbmi/vitessce/#readme).

## Install

1. Clone the repository

```sh
$ git clone git@github.com:haniffalab/sci-spatial-data.git
$ cd sci-spatial-data
```

2. Install nextflow by following the [official instruction](https://www.nextflow.io/index.html#GetStarted)
3. Install [Docker](https://docs.docker.com/engine/install/) and make sure it's in PATH
4. Install [Conda](https://docs.anaconda.com/anaconda/install/index.html) and make sure it's in PATH

Run
---

Create an [yaml file](template.yaml) for your datasets.

To convert h5ad to jsons:

```
nextflow run main.nf -params-file [your_params].yaml
```

To convert images to zarrs:

```
nextflow run main.nf -params-file [your_params].yaml -entry Image
nextflow run main.nf -params-file [your_params].yaml -entry Labels
```

To build a config file:

```
nextflow run main.nf -params-file [your_params].yaml -entry Config
```

Further reading:
--- 

Docker image pulling/local conda env creation are handled by nextflow. Please refer to [this](https://www.nextflow.io/docs/latest/getstarted.html) for detailed informtion.
