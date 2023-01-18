#! /bin/sh
#
# build-docker-imgs.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#

VERSION=0.0.1

docker build -t haniffalab/vitessce-pipeline-build-config:${VERSION} -f ./Dockerfile.build_config .
docker build -t haniffalab/vitessce-pipeline-processing:${VERSION} -f ./Dockerfile.processing .
docker build -t haniffalab/vitessce-pipeline-image-to-zarr:${VERSION} -f ./Dockerfile.image_to_zarr .
