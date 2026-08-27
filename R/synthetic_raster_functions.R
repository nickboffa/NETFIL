library(tidyverse)
library(terra)
library(gstat)
library(geodata)
library(sf)
library(geoR)

source("R/raster_pipeline_functions.R")   # allocate_pops, create_dist_csv, save_results, clean_raster

# Shared pipeline for fitting a Matern variogram to a real population raster
# and simulating synthetic ("SynthNNN") layouts from it. Used by
# R/automatic_variograms.R (the canonical Synth550 per country + diagnostics)
# and R/generate_synth550_realizations.R (additional named realizations for
# fit-sensitivity comparisons).

ISO_LOOKUP <- c(am_samoa = "ASM", samoa = "WSM", fiji = "FJI")

# 1998 (WSM) / most recent pre-model-period (ASM) population.
# https://data.worldbank.org/indicator/SP.POP.TOTL?locations=WS
TRUE_POP <- c(am_samoa = 54359, samoa = 178905)

fill_island_nas <- function(r, country, n_islands = 1) {
  iso <- ISO_LOOKUP[tolower(country)]
  if (is.na(iso)) stop("Unknown country '", country, "'. Options: ", paste(names(ISO_LOOKUP), collapse = ", "))

  boundary <- gadm(country = iso, level = 0, path = tempdir())
  parts    <- disagg(boundary)
  top_n    <- order(expanse(parts), decreasing = TRUE)[seq_len(n_islands)]
  islands  <- parts[top_n]

  mask_r   <- rasterize(islands, r, field = 1, background = NA)
  r_filled <- r
  r_filled[is.na(r_filled) & !is.na(mask_r)] <- 0
  r_filled
}

rescale_to_km <- function(r) {
  ext(r) <- ext(r) / 1000   # res() rescales automatically with a fixed cell count
  r
}

prep <- function(raster, fact = 1, zone = "EPSG:32702") {
  raster |>
    aggregate(fact = fact, fun = "sum", na.rm = TRUE) |>
    project(zone) |>
    rescale_to_km()
}

desspike_rank <- function(raster, tied_val = 0, radius = 0.2) {  # radius in km
  vals  <- values(raster, na.rm = FALSE)
  zeros <- which(vals == tied_val)

  focal_mean <- focal(raster, w = focalMat(raster, radius, "rectangle"),
                      fun = "mean", na.rm = TRUE)
  neigh_avg  <- values(focal_mean)[zeros]

  vals[zeros] <- tied_val + (rank(neigh_avg + rnorm(length(zeros), 0, 1e-6)) /
                               (length(zeros) + 1)) * 1e-6
  values(raster) <- vals
  raster
}

convert_normal_scores <- function(raster, method = "random") {
  vals   <- values(raster, na.rm = FALSE)
  not_na <- !is.na(vals)

  set.seed(1)
  ranked <- rank(vals[not_na], ties.method = method)
  ns     <- qnorm(ranked / (sum(not_na) + 1))

  ns_vals         <- vals
  ns_vals[not_na] <- ns
  values(raster)  <- ns_vals
  raster
}

fit_matern <- function(raster, cutoff = 45, width = 1.5, zone = 32702) {  # cutoff/width in km
  pts <- as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  colnames(pts) <- c("x", "y", "ns_value")
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = zone)

  v <- variogram(ns_value ~ 1, pts_sf, cutoff = cutoff, width = width)
  v_fit <- fit.variogram(v, vgm(psill = 1, model = "Mat"), fit.kappa = seq(0.01, 5, 0.01))
  list(v = v, v_fit = v_fit)
}

fit_country_variogram <- function(country_key, cfg) {
  prepped <- fill_island_nas(cfg$raster, country_key, n_islands = cfg$n_islands) |>
    prep(fact = cfg$prep_fact) |>
    desspike_rank(radius = 0.4) |>
    convert_normal_scores(method = "random")
  do.call(fit_matern, c(list(prepped), cfg$fit_args))
}

# Real per-cell population raster to sample synthetic layouts from. Scaled to
# the true country population here, before sampling — not after, on the
# small synthetic grid (which would inflate density: cramming the whole
# country's population into a synthetic footprint smaller than the real
# inhabited area). NA outside the island boundary (unsampleable); 0 for
# data-gap cells inside it (sampleable, same as any other land cell) —
# fill_island_nas() already draws that boundary, so we aggregate/scale
# without clean_raster()'s NA->0 step, which would otherwise turn ocean into
# sampleable zeros too.
real_population_raster <- function(country_key, cfg, sample_fact = 5, zone = "EPSG:32702") {
  agg <- fill_island_nas(cfg$raster, country_key, n_islands = cfg$n_islands) |>
    project(zone) |>
    aggregate(fact = sample_fact, fun = "sum", na.rm = TRUE)
  agg * (TRUE_POP[[country_key]] / global(agg, "sum", na.rm = TRUE)[1, 1])
}

# Sampling pool for make_synthetic_scale(): every land cell (0s included),
# NA (ocean) excluded.
real_cell_populations <- function(country_key, cfg, sample_fact = 5, zone = "EPSG:32702") {
  r <- real_population_raster(country_key, cfg, sample_fact = sample_fact, zone = zone)
  as.data.frame(r, na.rm = TRUE)[[1]]
}

# v_fit_model: a gstat vgm-class fit (e.g. fit_country_variogram(...)$v_fit),
# containing the Matern psill/range/kappa to simulate from.
make_synthetic_scale <- function(v_fit_model, n_side, real_vals, cell_km = 0.55,
                                  min_pop = 10, seed = 1) {
  mat      <- v_fit_model[v_fit_model$model == "Mat", ]
  cov_pars <- c(mat$psill, mat$range)

  set.seed(seed)
  fld <- grf(n = n_side^2, grid = "reg", nx = n_side, ny = n_side,
             xlims = c(0, (n_side - 1) * cell_km),
             ylims = c(0, (n_side - 1) * cell_km),
             cov.model = "matern", cov.pars = cov_pars, kappa = mat$kappa,
             messages = FALSE)

  # k-th largest real value -> k-th largest field cell
  n_grid  <- length(fld$data)
  sampled <- sample(real_vals, n_grid, replace = TRUE)
  ord     <- rank(fld$data, ties.method = "first")
  matched <- sort(sampled)[ord]

  r <- rast(data.frame(x = round(fld$coords[, 1], 6),
                       y = round(fld$coords[, 2], 6),
                       Population = matched),
            type = "xyz")

  # pop_target = grid's own current total: min_pop cleanup + integer rounding
  # only, no rescaling. Local origin — only relative positions matter downstream.
  set.seed(123)
  d  <- allocate_pops(r, min_pop = min_pop, pop_target = round(sum(matched)))
  x0 <- min(d$x); y0 <- min(d$y)
  d |>
    filter(Population >= min_pop) |>
    mutate(X = x - x0, Y = y - y0, Group = row_number()) |>
    select(Group, Population, X, Y)
}

find_n_side <- function(v_fit_model, target_groups, real_vals, start = 17,
                        cell_km = 0.55, min_pop = 10, seed = 1, max_side = 200, patience = 5) {
  best_n <- start
  best_dist <- Inf
  stall <- 0
  n <- start
  repeat {
    if (n > max_side) break
    dist_n <- abs(nrow(make_synthetic_scale(
      v_fit_model, n, real_vals, cell_km = cell_km, min_pop = min_pop, seed = seed
    )) - target_groups)
    if (dist_n < best_dist) {
      best_dist <- dist_n
      best_n <- n
      stall <- 0
    } else {
      stall <- stall + 1
      if (stall >= patience) break   # patience consecutive non-improving steps -> past the local minimum
    }
    n <- n + 1
  }
  best_n
}
