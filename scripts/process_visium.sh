#!/usr/bin/env bash
set -o errexit

. ./scripts/utils.sh

main() {
    # Download and process sample W3S whole embryo data. 
    # A single zarr store, and a cells file is created.

    get_CLI_args "$@"

    # Sample image data
    TIFF_IN="$INPUT/out.tif"
    ZARR_OUT="$OUTPUT/out.zarr"
    JSON_OUT="$OUTPUT/visium.raster.json"
    IMAGE_NAME="Visium Image Data"
}

### Main

main "$@"
