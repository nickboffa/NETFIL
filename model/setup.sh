#!/bin/bash

VALID_T2=(01 05 1 25)
VALID_SCALES=(Village Raster220 One Many Raster660)
T2=""; SCALE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t2)    T2="$2";    shift 2 ;;
        -scale) SCALE="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ " ${VALID_T2[*]} " =~ " $T2 " ]]        || { echo "Invalid -t2: $T2 (options: ${VALID_T2[*]})";        exit 1; }
[[ " ${VALID_SCALES[*]} " =~ " $SCALE " ]] || { echo "Invalid -scale: $SCALE (options: ${VALID_SCALES[*]})"; exit 1; }

echo "Cleaning and moving files"
./clean_inputs.sh

cp ../data/Scales/$SCALE/* ../data/
cp ../data/Fitted/$SCALE/Theta_$T2/TranParams.csv ../data/

echo "$SCALE" > ../data/current_scale.txt
