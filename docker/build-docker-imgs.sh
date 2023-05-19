#! /bin/sh
VERSION=0.3.1

docker build --platform=linux/amd64 -t haniffalab/webatlas-pipeline:${VERSION} -f ./Dockerfile .
