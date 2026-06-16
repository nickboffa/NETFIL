#!/bin/bash

declare -a array=(01 05 1 25)

if [[ ${array[*]} =~ $1 ]]
then
	echo Cleaning and moving files

	./clean_inputs.sh

	cp ../data/Scales/Raster220/groups.csv   ../data/groups.csv
	cp ../data/Scales/Raster220/euc_dist.csv ../data/euc_dist.csv
	cp ../data/Scales/Raster220/road_dist.csv ../data/road_dist.csv

	cp ../data/Fitted/Raster220/Theta_$1/Agg.txt    ../data/Fitted/
	cp ../data/Fitted/Raster220/Theta_$1/Theta1.txt  ../data/Fitted/
	cp ../data/Fitted/Raster220/Theta_$1/Work.txt    ../data/Fitted/
	cp ../data/Fitted/Raster220/Theta_$1/TranParams.csv ../data/
	echo "Raster220" > ../data/current_scale.txt
	exit 0
else
	echo "Incorrect input!"
	exit 1
fi
