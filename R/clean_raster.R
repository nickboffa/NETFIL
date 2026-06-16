library(terra)
library(tidyverse)

a110_raw <- rast("data/Rasters/am_samoa.tiff")

vil <- read_csv("data/villages_lat_lon.csv")
x_min <<- min(vil$X_Centroid)
y_min <<- min(vil$Y_Centroid)

Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj")

clean_raster <- function(raster, crs = "EPSG:32702", fact=1) {
  if (global(raster, fun="isNA") > 0) {
    raster <- subst(raster, NA, 0)
  }
  
  raster_utm <- project(raster, crs)
  
  raster_utm
}

create_raster_data <- function(raster, min_pop=2, pop_target = 54359) {

  data <- as.data.frame(raster, xy = TRUE, na.rm = FALSE) |>
    rename(Population = 3) # third column is band value
  
  data$Population[is.na(data$Population)] <- 0
  data$Population <- (pop_target / sum(data$Population)) * data$Population
  
  data |> 
    filter(Population >= min_pop) |> 
    mutate(
      Population = round(Population),
      X = x - x_min,
      Y = y - y_min,
      Group = row_number()
    ) |>
    select(Group, Population, X, Y)
}

create_euc_dist <- function(groups) {
  coords <- as.matrix(groups[, c("X", "Y")])
  d <- round(as.matrix(dist(coords)))
  euc_dists <- as.data.frame(d) |>
    mutate(
      X = groups$Group
    ) |>
    relocate(X)
}

create_road_dist <- function(groups) {
  coords <- as.matrix(groups[, c("X", "Y")])
  d <- round(as.matrix(dist(coords, method = "manhattan")))
  road_dists <- as.data.frame(d) |>
    mutate(
      X = groups$Group
    ) |>
    relocate(X)
}

a110 <- clean_raster(a110_raw)
data_a110 <- create_raster_data(a110, pop_target=52380)
euc_a110 <- create_euc_dist(data_a110)
road_a110 <- create_road_dist(data_a110)

a220 <- aggregate(a110, fact=2, fun="sum")
data_a220 <- create_raster_data(a220, pop_target=52380)
euc_a220 <- create_euc_dist(data_a220)
road_a220 <- create_road_dist(data_a220)

a330 <- aggregate(a110, fact=3, fun="sum")
data_a330 <- create_raster_data(a330, pop_target=52380)
euc_a330 <- create_euc_dist(data_a330)

a440 <- aggregate(a110, fact=4, fun="sum")
data_a440 <- create_raster_data(a440, pop_target=52380)
euc_a440 <- create_euc_dist(data_a440)

a550 <- aggregate(a110, fact=5, fun="sum")
data_a550 <- create_raster_data(a550, pop_target=52380)
euc_a550 <- create_euc_dist(data_a550)

a660 <- aggregate(a110, fact=6, fun="sum")
data_a660 <- create_raster_data(a660, pop_target=52380)
euc_a660 <- create_euc_dist(data_a660)
road_a660 <- create_road_dist(data_a660)

write_csv(data_a220, "data/Scales/Raster220/groups.csv")
write_csv(euc_a220, "data/Scales/Raster220/euc_dist.csv")
write_csv(road_a220, "data/Scales/Raster220/road_dist.csv")

write_csv(data_a660, "data/Scales/Raster660/groups.csv")
write_csv(euc_a660, "data/Scales/Raster660/euc_dist.csv")
write_csv(road_a660, "data/Scales/Raster660/road_dist.csv")


global(a110, fun=\(x) sum(x >= 1, na.rm=TRUE))

global(a660, fun=\(x) sum(x >= 1, na.rm=TRUE))
length(read.csv("data/Scales/Many/groups.csv")$Population)
length(read.csv("data/Scales/One/groups.csv")$Population)
length(read.csv("data/Scales/Village/groups.csv")$Population)


sum(read.csv("data/Scales/Many/groups.csv")$Population)
sum(read.csv("data/Scales/One/groups.csv")$Population)
sum(read.csv("data/Scales/Village/groups.csv")$Population)

sum(data_a660$Population)
