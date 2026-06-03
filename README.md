# NETFIL

This version edited from Callum's.

## Coordinates

The coordinates are UTM (presumably the AS zone), subtracting the minimum village coordinates. This means (X, Y) = (0,0) is the south-western corner. The absolute UTM easting is 518002.3, and a northing of 8412742.7.

## Data

Distances are in metres (e.g. in euc_dist.csv)

## Data

American Samoa population density data downloaded from

https://geonode.pacificdata.org/catalogue/#/dataset/1027
https://pacificdata.org/data/dataset/asm-pop-grid-2020-1027

## Mapping paper to code

paper | code
theta1

k | agg_param

w | worktonot

## Issues

Population in One and Village is 54359

But population in Many is 52380.

For the moment I have made Raster scale to have the population of 54359, and then rounded. This made

> sum(data_a220$Population)
[1] 54180
> sum(data_a110$Population)
[1] 53482


Is it better to filter out small population cells, or should I really combine them into larger ones?

Currently commuting prop is fixed at 50%
As cell size increases/decreases, this should change

Could make commuting prop a fitted parameter, fitting to ICC (can calculate ICC of a simulation by combining groups to form a village-sized group, then computing ICC. Can therefore make it a part of ABC.

Similarly, should the distribution of prevalences among groups also change?

Antigen positive (cout-ed when running model) is a double when it should be an int

Running the model should make a .netfil file (like a .iqtree file) that contains the command run, the start and end time, a copy of MDA_params.csv and InitParams.csv
Maybe params.h?
Running model will create a folder within output/ not a single csv

Values of theta1, theta2, agg, worktonot used
Etc.

sim_years should depend on ABC_fitting = true or note (7 if true, 21 otherwise)

Check whether the population does grow or shrink – it seems that it doesn’t?

What to do with single-person groups? And with rounding?


Currently:
    int recalc_years = 100; //how often we want to recalc commuters
    char distance_type = 'r'; // r for road distance, e for euclidean 

But doesn’t the paper say it should use euclidean distance? And regardless, this should be in param.h surely?


In ABC, am unsuer what epsilon was used, and the parameter selected seems to be slightly different to the median? It's not the mean or mode though.