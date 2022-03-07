#!/usr/bin/env bash
set -o errexit

. ./scripts/utils.sh

main() {
    # Download and process sample W3S whole embryo data. 
    # A single zarr store, and a cells file is created.

    get_CLI_args "$@"

    # Single cell data
    H5AD_IN="$INPUT/data.h5ad"
    CELLS_OUT="$OUTPUT/cells.json"
    CELL_SETS_OUT="$OUTPUT/cell-sets.json"
    MATRIX_OUT="$OUTPUT/clusters.json"

    if [ -e "$CELLS_OUT" ]
    then
        echo "Skipping cells -- output already exists: $CELLS_OUT"
    else
        echo 'Generating cells JSON may take a while...'
        CMD="$BASE/python/spot_reader.py
            --h5ad_file $H5AD_IN
            --cells_file $CELLS_OUT
            --cell_sets_file $CELL_SETS_OUT
            --matrix_file $MATRIX_OUT"
        echo "Running: $CMD"
        eval $CMD
    fi   
}

### Main

main "$@"
