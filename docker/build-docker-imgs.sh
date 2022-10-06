#! /bin/sh
#
# build-docker-imgs.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#

VERSION=0.0.1

docker build -t hamat/webatlas-zarr:${VERSION} -f ./Dockerfile.zarr .
docker build -t hamat/webatlas-build-config:${VERSION} -f ./Dockerfile.build_config .
docker build -t hamat/webatlas-generate-label:${VERSION} -f ./Dockerfile.generate_label .
docker build -t hamat/webatlas-ome-zarr-metadata:${VERSION} -f ./Dockerfile.ome_zarr_metadata .
docker build -t hamat/webatlas-router:${VERSION} -f ./Dockerfile.router .
docker build -t hamat/webatlas-image-to-zarr:${VERSION} -f ./Dockerfile.image_to_zarr .
