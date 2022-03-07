#!/usr/bin/env bash
set -o errexit

. ./scripts/utils.sh

main() {
    # Download and process iss data. The following outputs are produced: 
    # - Microscopy zarr store
    # - Cell boundary zarr store
    # - cell.json (http://vitessce.io/docs/data-file-types/#cellsjson)
    # - cell-sets.json (http://vitessce.io/docs/data-file-types/#cell-setsjson)
    # - raster.json (http://vitessce.io/docs/data-file-types/#rasterjson)
    # - clusters.json (http://vitessce.io/docs/data-file-types/#clustersjson)
    # - genes.json (http://vitessce.io/docs/data-file-types/#genesjson)

    get_CLI_args "$@"

    # Image data
    TIFF_IN="$INPUT/imaging-data.tif"
    ZARR_OUT="$OUTPUT/out.zarr"
    JSON_OUT="$OUTPUT/raster.json"
    IMAGE_NAME="Image Data"

    # Single cell data
    H5AD_IN="$INPUT/sc-data.h5ad"
    CELLS_OUT="$OUTPUT/cells.json"
    CELL_SETS_OUT="$OUTPUT/cell-sets.json"
    MATRIX_OUT="$OUTPUT/clusters.json"

    # Cloud source
    SOURCE_URL="https://storage.googleapis.com/webatlas-embryo/sample"
    SOURCE_BUCKET="webatlas-embryo"

    # Cloud target
    RELEASE=${CLOUD_TARGET//webatlas-vitessce-data\//}
    DEST_URL="https://storage.googleapis.com/$RELEASE/sample/"

    echo "Download and process sample image data..."

    # [ -e "$H5AD_IN" ] || \
    #     #gsutil cp gs://"$SOURCE_BUCKET"/sample/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys.h5ad "$H5AD_IN"
    #     #wget "$SOURCE_URL/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys.h5ad" -O "$H5AD_IN"
    #     #wget "https://a04fcc815aa920b9c7e028bb79f7c2db29d0682c.cog.sanger.ac.uk/39bfafda600ff69122887bce04f4efb88f767caa/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys_shift_20_$i.h5ad" -O "$H5AD_IN"

    if [ -e "$ZARR_OUT" ]
    then
        echo "Skipping tiling -- output already exists: $ZARR_OUT"
    else
        if [ -e "$TIFF_IN" ]
        then
          echo "Not copying $TIFF_IN from s3 - already exists or testing"
        # else
        #     gsutil cp gs://"$SOURCE_BUCKET"/sample/out.tif "$TIFF_IN"
        #     wget "$SOURCE_URL/out.tif" -O "$TIFF_IN"
        fi
        echo 'Converting OME-TIFF to zarr may take a while...'
        CMD="$BASE/python/ome_tiff_reader.py
            --input_tiff $TIFF_IN
            --output_zarr $ZARR_OUT
            --image_json $JSON_OUT
            --image_name '$IMAGE_NAME'
            --dest_url $DEST_URL"
        echo "Running: $CMD"
        eval $CMD
    fi

    echo "Download image data..."

    # [ -e "$H5AD_IN" ] || \
    #     #gsutil cp gs://"$SOURCE_BUCKET"/sample/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys.h5ad "$H5AD_IN"
    #     #wget "$SOURCE_URL/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys.h5ad" -O "$H5AD_IN"
    #     #wget "https://a04fcc815aa920b9c7e028bb79f7c2db29d0682c.cog.sanger.ac.uk/39bfafda600ff69122887bce04f4efb88f767caa/N1234F_OB10037_x20_GMMdecoding_shift_sep3_cutoff-11_doublet_with_ploys_shift_20_$i.h5ad" -O "$H5AD_IN"

    if [ -e "$CELLS_OUT" ]
    then
        echo "Skipping cells -- output already exists: $CELLS_OUT"
    else
        echo 'Generating cells JSON may take a while...'
        CMD="$BASE/python/h5ad_reader.py
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
