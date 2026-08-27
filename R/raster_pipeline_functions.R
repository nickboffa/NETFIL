# Shared functions for turning a population raster (real or synthetic) into
# the groups.csv / euc_dist / road_dist files a scale needs. Used by both
# R/create_required_raster_csvs.R (real rasters) and R/automatic_variograms.R
# (synthetic rasters) — source() this file rather than duplicating these.

library(terra)
library(tidyverse)

clean_raster <- function(raster, crs = "EPSG:32702") {
  if (global(raster, fun = "isNA") > 0) {
    raster <- subst(raster, NA, 0)
  }
  project(raster, crs)
}

# Rescale a population raster to pop_target, redistributing cells that fall
# below min_pop geographically: each excluded cell donates 100% of its
# population to the nearest still-included cell (Euclidean distance).
# Tie-break 1: lowest current population. Tie-break 2: random.
# Returns a data frame with columns x, y, Population (all cells; excluded → 0).
allocate_pops <- function(raster, min_pop = 10, pop_target) {
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

create_raster_data <- function(raster, min_pop = 10, pop_target, x_min, y_min) {
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
groups_to_raster <- function(groups, template, value = "Population", x_offset, y_offset) {
  pts <- groups |>
    mutate(x = X + x_offset, y = Y + y_offset) |>
    rename(val = !!sym(value))

  v <- vect(pts, geom = c("x", "y"), crs = crs(template))
  rasterize(v, template, field = "val", background = NA)
}

# "manhattan" = L1/road, "euclidean" = L2/'euc'
create_dist_csv <- function(groups, dist_method = "euclidean") {
  coords <- as.matrix(groups[, c("X", "Y")])
  d <- round(as.matrix(dist(coords, method = dist_method)))

  as.data.frame(d) |>
    mutate(X = groups$Group) |>
    relocate(X)
}

write_dist_binary <- function(dist_df, path) {
  mat  <- as.matrix(dist_df[, -1])      # drop group-ID column; n×n numeric matrix
  vals <- mat[lower.tri(mat)]            # lower-tri col-major = upper-tri row-major (symmetric)
  con  <- file(path, "wb")
  writeBin(as.double(vals), con)
  close(con)
}

save_results <- function(res, country = c("WSM", "ASM"), name) {
  base <- file.path("data", country, "Scales", name)
  dir.create(base, showWarnings = FALSE, recursive = TRUE)
  write_csv(res$data, file.path(base, "groups.csv"))
  write_csv(res$euc,  file.path(base, "euc_dist.csv"))
  write_csv(res$road, file.path(base, "road_dist.csv"))
  write_dist_binary(res$euc,  file.path(base, "euc_dist.bin"))
  write_dist_binary(res$road, file.path(base, "road_dist.bin"))
}
