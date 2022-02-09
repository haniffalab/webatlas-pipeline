#!/usr/bin/env bash

die() { set +v; echo "$*" 1>&2 ; exit 1; }

add_CLI_ARGS() {
    # Helper for process_cells to build argument list.

    FILE_TYPE=$1
    DATA_TITLE=$2

    FILE="$OUTPUT/$DATA_TITLE.$FILE_TYPE.json"
    if [ -e "$FILE" ]
    then
        echo "$FILE_TYPE output already exists: $FILE"
    else
        CLI_ARGS="$CLI_ARGS --${FILE_TYPE}_file $FILE"
    fi
}

usage() {
    [ -z "$1" ] || echo "$1 is not a directory." 1>&2
    die "Usage: $0 -b <directory> -i <directory> -o <directory> -t <target>
    -b   Base directory
    -i   Input directory
    -o   Output directory
    -t   Cloud target"
}

get_CLI_args(){
    # echo "Parsing: $@"
    while getopts "b:i:o:t:" arg; do
        count=$(($count + 1))
        case $arg in
            b)
                [ -d "$OPTARG" ] || usage "$OPTARG"
                BASE=${OPTARG}
                ;;
            i)
                [ -d "$OPTARG" ] || usage "$OPTARG"
                INPUT=${OPTARG}
                ;;
            o)
                [ -d "$OPTARG" ] || usage "$OPTARG"
                OUTPUT=${OPTARG}
                ;;
            t)
                CLOUD_TARGET=${OPTARG}
                ;;
        esac
    done

    if [ "$count" != "4" ]
    then
        echo $count
        echo "All 4 arguments are required"
        usage
    fi
}
