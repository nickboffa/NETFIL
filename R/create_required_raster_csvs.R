library(terra)
library(tidyverse)
select <- dplyr::select

a110_raw <- rast("data/Rasters/ASM.tiff")
s110_raw <- rast("data/Rasters/WSM.tif")

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

# Rescale a population raster to pop_target, redistributing cells that fall
# below min_pop geographically: each excluded cell donates 100% of its
# population to the nearest still-included cell (Euclidean distance).
# Tie-break 1: lowest current population. Tie-break 2: random.
# Returns a data frame with columns x, y, Population (all cells; excluded → 0).
allocate_pops <- function(raster, min_pop = 10, pop_target = 54359) {
  data <- as.data.frame(raster, xy = TRUE, na.rm = FALSE)
  names(data)[3] <- "Population"
  data$Population[is.na(data$Population)] <- 0

  pops   <- data$Population
  coords <- as.matrix(data[, c("x", "y")])

  pos       <- pops > 0
  exact_pop <- numeric(length(pops))
  exact_pop[pos] <- (pop_target / sum(pops[pos])) * pops[pos]

  included <- exact_pop >= min_pop
  to_give  <- which(exact_pop > 0 & !included)

  for (i in to_give) {
    incl_idx   <- which(included)
    dx         <- coords[incl_idx, 1] - coords[i, 1]
    dy         <- coords[incl_idx, 2] - coords[i, 2]
    dists      <- sqrt(dx^2 + dy^2)
    min_d      <- min(dists)
    candidates <- incl_idx[dists == min_d]

    if (length(candidates) > 1) {
      min_p      <- min(exact_pop[candidates])
      candidates <- candidates[exact_pop[candidates] == min_p]
    }
    if (length(candidates) > 1) candidates <- sample(candidates, 1)

    exact_pop[candidates] <- exact_pop[candidates] + exact_pop[i]
    exact_pop[i] <- 0
  }

  # Integer allocation: floor + distribute remainders by largest fractional part
  valid         <- exact_pop > 0
  floor_pop     <- floor(exact_pop)
  rem_rank      <- rep(Inf, length(exact_pop))
  rem_rank[valid] <- rank(-(exact_pop[valid] - floor_pop[valid]), ties.method = "random")
  n_to_ceil     <- pop_target - sum(floor_pop)

  data$Population <- floor_pop + (rem_rank <= n_to_ceil)
  data
}

create_raster_data <- function(raster, min_pop = 10, pop_target = 54359) {
  set.seed(123)
  data <- allocate_pops(raster, min_pop = min_pop, pop_target = pop_target)

  data |>
    filter(Population >= min_pop) |>
    mutate(
      X = x - x_min,
      Y = y - y_min,
      Group = row_number()
    ) |>
    select(Group, Population, X, Y)
}

# Reconstruct a plottable SpatRaster from groups data.
# `template` should be the cleaned raster (or an aggregate of it) that was
# passed to create_raster_data — it supplies the resolution, extent, and CRS.
# `value` is the column in `groups` to use as the raster value (default: Population).
groups_to_raster <- function(groups, template,
                              value = "Population",
                              x_offset = x_min, y_offset = y_min) {
  pts <- groups |>
    mutate(x = X + x_offset, y = Y + y_offset) |>
    rename(val = !!sym(value))

  v <- vect(pts, geom = c("x", "y"), crs = crs(template))
  rasterize(v, template, field = "val", background = NA)
}

# "manhattan" = L1/road
# "euclidean" = L2/'euc'

create_dist_csv <- function(groups, dist_method = "euclidean") {
  coords <- as.matrix(groups[, c("X", "Y")])
  d <- round(as.matrix(dist(coords, method=dist_method)))
  
  as.data.frame(d) |>
    mutate(
      X = groups$Group
    ) |>
    relocate(X)
}

raster <- a110_raw

facts <- c(1, 2, 3, 4, 5, 6)

rasters <- list(raster |> clean_raster())
for (i in 2:max(facts)) {
  rasters[[i]] <- aggregate(rasters[[1]], fact=facts[i], fun="sum")
}

# samoan pop target from
# https://data.worldbank.org/indicator/SP.POP.TOTL?locations=WS
results_orig <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 10)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

lapply(seq_along(rasters), function(i) {
  cat(i, nrow(results_orig[[i]]$data), sum(results_orig[[i]]$data$Population), "\n")
})


save_results <- function(res, country = c("WSM", "ASM"), name) {
  base <- file.path("data", country, "Scales", name)
  write_csv(res$data, file.path(base, "groups.csv"))
  write_csv(res$euc, file.path(base, "euc_dist.csv"))
  write_csv(res$road, file.path(base, "road_dist.csv"))
}




save_results(results_orig[[2]], "Raster220")

save_results(results_orig[[5]], "Raster550")

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
  data <- create_raster_data(rasters[[i]], min_pop = 2, pop_target = 10000)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

save_results(results_10k[[2]], "Raster220_10k")
save_results(results_10k[[6]], "Raster660_10k")

results_20k <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 4, pop_target = 20000)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

save_results(results_20k[[2]], "Raster220_20k")
save_results(results_20k[[6]], "Raster660_20k")

results_100k <- lapply(seq_along(rasters), function(i) {
  data <- create_raster_data(rasters[[i]], min_pop = 4, pop_target = 100000)
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

s_results <- lapply(seq_along(s_rasters[2:9]), function(i) {
  data <- create_raster_data(s_rasters[[i+1]], min_pop = 90, pop_target = 193023)
  list(
    data = data,
    euc  = create_dist_csv(data, "euclidean"),
    road = create_dist_csv(data, "manhattan")
  )
})

for (i in seq_along(s_rasters)) {
  cat(i+1, nrow(s_results[[i]]$data), sum(s_results[[i]]$data$Population), "\n")
}

data_wsm <- create_raster_data(rasters[[5]], min_pop = 30, pop_target = 178905)

r_wsm <- groups_to_raster(data_wsm, rasters[[5]])


subst(rasters[[5]], 0, NA) |> plot()
plot(r_wsm)

res_wsm <- list(
  data = data_wsm,
  euc  = create_dist_csv(data_wsm, "euclidean"),
  road = create_dist_csv(data_wsm, "manhattan")
)

