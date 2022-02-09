#!/usr/bin/env bash
set -o errexit

die() { set +v; echo "$*" 1>&2 ; exit 1; }

if [ "$#" -ne 0 ]; then
    die 'Collects source data, processes it, and pushes to S3. No commandline arguments, but looks for two environment variables:
  - "NO_PUSH=true" will fetch and process data, but not push to S3.
  - "CI=true" will use fixtures, rather than fetching, and also will not push to S3.'
fi

main() {
    INPUT="$FILES/input"
    OUTPUT="$FILES/output"

    for DATASET in iss visium; do
        INPUT_SET="$INPUT/$DATASET"
        OUTPUT_SET="$OUTPUT/$DATASET"
        [ -d "$INPUT_SET" ] || mkdir -p "$INPUT_SET"
        [ -d "$OUTPUT_SET" ] || mkdir -p "$OUTPUT_SET"

        echo "Processing $DATASET ..."
        ./scripts/process_$DATASET.sh \
            -b "$BASE" \
            -i "$INPUT_SET" \
            -o "$OUTPUT_SET" \
            -t "$CLOUD_TARGET"
    done
}

### Globals

BASE=`pwd`
CLOUD_TARGET=`cat cloud_target.txt`
FILES="$BASE/data-files"

### Main

main "$FILES"
