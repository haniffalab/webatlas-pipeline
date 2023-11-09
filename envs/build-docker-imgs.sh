#! /bin/sh
VERSION=0.4.1

docker build --platform=linux/amd64 -t haniffalab/webatlas-pipeline:${VERSION} -f ./Dockerfile .
cd build_config/
docker build --platform=linux/amd64 -t haniffalab/webatlas-pipeline-build_config:${VERSION} -f ./Dockerfile .
