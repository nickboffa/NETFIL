library(terra)
library(tidyverse)
select <- dplyr::select

source("R/raster_pipeline_functions.R")   # clean_raster, allocate_pops, create_raster_data,
                                           # groups_to_raster, create_dist_csv,
                                           # write_dist_binary, save_results

a110_raw <- rast("data_storage/Rasters/ASM.tiff")
s110_raw <- rast("data_storage/Rasters/WSM.tif")

vil <- read_csv("data_storage/villages_lat_lon.csv")
x_min <- min(vil$X_Centroid)
y_min <- min(vil$Y_Centroid)

Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj")



raster <- a110_raw

facts <- c(1, 2, 3, 4, 5, 6)

rasters <- list(raster |> clean_raster())
for (i in 2:max(facts)) {
  rasters[[i]] <- aggregate(rasters[[1]], fact=facts[i], fun="sum")
}

# samoan pop target from
# https://data.worldbank.org/indicator/SP.POP.TOTL?locations=WS
results_orig <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 10, pop_target = 54359, x_min = x_min, y_min = y_min)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

lapply(seq_along(rasters), function(i) {
  cat(i, nrow(results_orig[[i]]$data), sum(results_orig[[i]]$data$Population), "\n")
})

save_results(results_orig[[2]], "ASM", "Raster220")

save_results(results_orig[[5]], "ASM", "Raster550")

save_results(results_orig[[6]], "ASM", "Raster660")

save_results(res_wsm, "WSM", "Raster550")

# Previously

# > sum(read.csv("data/Scales/Many/groups.csv")$Population)
# [1] 52380
# > sum(read.csv("data/Scales/One/groups.csv")$Population)
# [1] 54359
# > sum(read.csv("data/Scales/Village/groups.csv")$Population)
# [1] 54359

# Note <-  Many was proportionally rescaled, with min_pop 0
# to get to the same population as the others.
many <- read.csv("data/Scales/Many/groups.csv")
many$Population <- rescale_pops(many$Population, min_pop = 0)
write_csv(many, "data/Scales/Many/groups.csv")

length(read.csv("data/Scales/Many/groups.csv")$Population)
length(read.csv("data/Scales/One/groups.csv")$Population)
length(read.csv("data/Scales/Village/groups.csv")$Population)

sum(read.csv("data/Scales/Many/groups.csv")$Population)
sum(read.csv("data/Scales/One/groups.csv")$Population)
sum(read.csv("data/Scales/Village/groups.csv")$Population)


results_10k <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 2, pop_target = 10000, x_min = x_min, y_min = y_min)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

save_results(results_10k[[2]], "Raster220_10k")
save_results(results_10k[[6]], "Raster660_10k")

results_20k <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 4, pop_target = 20000, x_min = x_min, y_min = y_min)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

save_results(results_20k[[2]], "Raster220_20k")
save_results(results_20k[[6]], "Raster660_20k")

results_100k <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 4, pop_target = 100000, x_min = x_min, y_min = y_min)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

save_results(results_100k[[2]], "Raster220_100k")
save_results(results_100k[[6]], "Raster660_100k")


for (i in seq_along(rasters)) {
  cat(i, nrow(results_orig[[i]]$data), sum(results_orig[[i]]$data$Population), "\n")
}

for (i in seq_along(rasters)) {
  cat(i, nrow(results_10k[[i]]$data), sum(results_10k[[i]]$data$Population), "\n")
}






hist110 <- data_a110 |> 
  ggplot(aes(Population)) +
  geom_histogram(binwidth=1, center = 0,
                 fill="grey", color="black")

results[[2]]$data |> 
  ggplot(aes(Population)) +
  geom_histogram(binwidth=1, center = 0,
                 fill="grey", color="black")

results[[4]]$data |> 
  ggplot(aes(Population)) +
  geom_histogram(binwidth=1, center = 0,
                 fill="grey", color="black")

## SAMOA

s_rasters <- list(s110_raw |> clean_raster())
for (i in 2:5) {
  s_rasters[[i]] <- aggregate(s_rasters[[1]], fact=i, fun="sum")
}

data_wsm <- create_raster_data(s_rasters[[5]], min_pop = 10, pop_target = 178905, x_min = x_min, y_min = y_min)

r_wsm <- groups_to_raster(data_wsm, s_rasters[[5]], x_offset = x_min, y_offset = y_min)


subst(s_rasters[[5]], 0, NA) |> plot()
plot(r_wsm)

res_wsm <- list(
  data = data_wsm,
  euc  = create_dist_csv(data_wsm, "euclidean"),
  road = create_dist_csv(data_wsm, "manhattan")
)

save_results(res_wsm, "WSM", "Raster550")

