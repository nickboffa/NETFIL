#!/bin/bash

VALID_T2=(01 05 1 25)
VALID_SCALES=(Village Raster220 One Many Raster660 Raster550)
VALID_COUNTRIES=(ASM WSM)
T2=""; SCALE=""; COUNTRY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t2)      T2="$2";      shift 2 ;;
        -scale)   SCALE="$2";   shift 2 ;;
        -country) COUNTRY="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ " ${VALID_T2[*]} " =~ " $T2 " ]]           || { echo "Invalid -t2: $T2 (options: ${VALID_T2[*]})";           exit 1; }
[[ " ${VALID_SCALES[*]} " =~ " $SCALE " ]]    || { echo "Invalid -scale: $SCALE (options: ${VALID_SCALES[*]})";    exit 1; }
[[ " ${VALID_COUNTRIES[*]} " =~ " $COUNTRY " ]] || { echo "Invalid -country: $COUNTRY (options: ${VALID_COUNTRIES[*]})"; exit 1; }

echo "Cleaning and moving files"
./clean_inputs.sh

cp ../data/$COUNTRY/Scales/$SCALE/* ../data/
cp ../data/$COUNTRY/Fitted/$SCALE/Theta_$T2/tran_params.csv ../data/
cp ../data/$COUNTRY/mda_params.csv ../data/

echo "$SCALE" > ../data/current_scale.txt
echo "$COUNTRY" > ../data/current_country.txt
