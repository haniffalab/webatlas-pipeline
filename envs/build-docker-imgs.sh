#! /bin/sh
VERSION=0.5.0

docker build --platform=linux/amd64 -t haniffalab/webatlas-pipeline:${VERSION} -f ./Dockerfile .
cd build_config/
docker build --platform=linux/amd64 -t haniffalab/webatlas-pipeline-build-config:${VERSION} -f ./Dockerfile .
