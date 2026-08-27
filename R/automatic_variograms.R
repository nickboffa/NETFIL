Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")
library(patchwork)

source("R/synthetic_raster_functions.R")

# https://geonode.pacificdata.org/catalogue/#/dataset/1027
# https://pacificdata.org/data/dataset/asm-pop-grid-2020-1027
am_samoa <- rast("data_storage/Rasters/ASM.tiff")

# https://pacificdata.org/data/dataset/population-grid-samoa-2025
samoa <- rast("data_storage/Rasters/WSM.tif")

# Cuts a sliver of stray raster cells far from Tutuila that confuses the
# single-island mask below.
x_lon <- init(am_samoa, "x")
x_lat <- init(am_samoa, "y")
am_samoa[x_lon > -170.58 & x_lat < -14.275] <- NA

countries <- list(
  am_samoa = list(raster = am_samoa, label = "American Samoa", n_islands = 1,
                  prep_fact = 2, fit_args = list(cutoff = 20, width = 0.5)),
  samoa    = list(raster = samoa, label = "Samoa", n_islands = 2,
                  prep_fact = 2, fit_args = list())
)

results <- lapply(names(countries), function(name) {
  list(v_fit = fit_country_variogram(name, countries[[name]]), label = countries[[name]]$label)
})
names(results) <- names(countries)

format_vfit_title <- function(country_name, v_fit_obj) {
  fit <- v_fit_obj$v_fit
  mat_row <- fit[fit$model != "Nug", ]
  nug_row <- fit[fit$model == "Nug", ]

  nugget_str <- if (nrow(nug_row) > 0 && nug_row$psill > 0)
    sprintf(", nugget=%.3f", nug_row$psill) else ""

  subtitle   <- sprintf("psill=%.3f, range=%.2fkm, kappa=%.2f%s",
                        mat_row$psill, mat_row$range, mat_row$kappa, nugget_str)
  show(plot(v_fit_obj$v, v_fit_obj$v_fit, main = sprintf("%s\n%s", country_name, subtitle)))
}

for (name in names(results)) {
  format_vfit_title(results[[name]]$label, results[[name]]$v_fit)
}

# --------------------------------------------------------------------------
# Diagnostic: how flat is the kappa fit? SSE (fit.variogram()'s WLS
# objective) at each candidate kappa, fixing kappa and refitting psill+range
# each time -- the same grid fit_matern() already does internally, just with
# every candidate's SSE exposed instead of only the argmin. A flat curve
# near the minimum means kappa isn't well pinned down by the data.
# --------------------------------------------------------------------------

kappa_sse_curve <- function(v, kappa_seq = seq(0.01, 5, 0.01)) {
  map_dfr(kappa_seq, function(k) {
    fit <- tryCatch(
      suppressWarnings(fit.variogram(v, vgm(psill = 1, model = "Mat", kappa = k), fit.kappa = FALSE)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    tibble(kappa = k, sse = attr(fit, "SSErr"))
  })
}

for (name in names(results)) {
  v          <- results[[name]]$v_fit$v
  mat_row    <- results[[name]]$v_fit$v_fit
  chosen_k   <- mat_row$kappa[mat_row$model == "Mat"]
  curve_df   <- kappa_sse_curve(v)

  p <- ggplot(curve_df, aes(kappa, sse)) +
    geom_line() +
    geom_vline(xintercept = chosen_k, linetype = "dashed", color = "red") +
    labs(title = sprintf("%s: SSE across kappa grid search", results[[name]]$label),
         subtitle = sprintf("chosen kappa = %.2f", chosen_k),
         x = "kappa", y = "SSE (gstat fit.variogram, fixed kappa)") +
    theme_minimal()
  show(p)
}

# --------------------------------------------------------------------------
# Build + save Synth550 for every country, all targeting the same group
# count (ASM's real Raster550 count) rather than each country's own real
# count — the point of Synth550 for WSM is a smaller, faster-to-run
# population, not a size-matched replica of WSM's real (1694-group) Raster550.
# --------------------------------------------------------------------------

TARGET_GROUPS <- nrow(read_csv("data/ASM/Scales/Raster550/groups.csv", show_col_types = FALSE))

synth_results <- lapply(names(countries), function(key) {
  iso <- ISO_LOOKUP[[key]]
  cfg <- countries[[key]]
  v_fit_model <- results[[key]]$v_fit$v_fit
  real_vals   <- real_cell_populations(key, cfg)

  message(sprintf("%s: sweeping n_side to hit ~%d groups...", iso, TARGET_GROUPS))
  n_best <- find_n_side(v_fit_model, TARGET_GROUPS, real_vals)
  groups <- make_synthetic_scale(v_fit_model, n_best, real_vals)
  message(sprintf("  n_side=%d -> %d groups (target %d)", n_best, nrow(groups), TARGET_GROUPS))

  res <- list(
    data = groups,
    euc  = create_dist_csv(groups, "euclidean"),
    road = create_dist_csv(groups, "manhattan")
  )
  save_results(res, iso, "Synth550")
  res
})
names(synth_results) <- names(countries)

# --------------------------------------------------------------------------
# Additional named realizations (Synth550_v1, _v2, ...), distinct from the
# canonical "Synth550" above (seed=1) — for measuring how much ABC-fit
# parameters/predictions depend on which random synthetic layout was drawn,
# not on stochastic disease dynamics. Reuses the variogram fit and real_vals
# pool already computed above, so no re-fitting cost per extra realization.
# clustering_params.csv is copied from the country's real Raster550, same
# reasoning as the canonical version's transmission params being copied —
# a synthetic layout has no real ICC data of its own to derive one from.
# --------------------------------------------------------------------------

EXTRA_REALIZATIONS <- list(samoa = 1)   # country_key -> how many extra versions to build

for (key in names(EXTRA_REALIZATIONS)) {
  iso <- ISO_LOOKUP[[key]]
  v_fit_model <- results[[key]]$v_fit$v_fit
  real_vals   <- real_cell_populations(key, countries[[key]])
  real_clustering_csv <- file.path("data", iso, "Scales", "Raster550", "clustering_params.csv")

  for (i in seq_len(EXTRA_REALIZATIONS[[key]])) {
    seed <- i + 1   # +1 so v1 doesn't collide with the canonical Synth550's seed=1
    name <- sprintf("Synth550_v%d", i)

    message(sprintf("%s: sweeping n_side for %s (seed=%d)...", iso, name, seed))
    n_best <- find_n_side(v_fit_model, TARGET_GROUPS, real_vals, seed = seed)
    groups <- make_synthetic_scale(v_fit_model, n_best, real_vals, seed = seed)
    message(sprintf("  n_side=%d -> %d groups (target %d)", n_best, nrow(groups), TARGET_GROUPS))

    res <- list(
      data = groups,
      euc  = create_dist_csv(groups, "euclidean"),
      road = create_dist_csv(groups, "manhattan")
    )
    save_results(res, iso, name)

    if (file.exists(real_clustering_csv)) {
      file.copy(real_clustering_csv, file.path("data", iso, "Scales", name, "clustering_params.csv"),
                overwrite = TRUE)
    }
  }
}

# --------------------------------------------------------------------------
# Diagnostic: real vs synthetic layout, per country
# --------------------------------------------------------------------------

to_km <- function(df, cell_m) {
  df |> mutate(X_km = (X - min(X)) / 1000, Y_km = (Y - min(Y)) / 1000)   # X/Y in metres -> km
}

# Both panels get the same xlim/ylim/fill limits, so a physically smaller or
# sparser synthetic layout actually looks smaller/sparser on the page —
# unlike col/row grid-index axes, which auto-scale each panel to fill the
# same plot area regardless of real physical extent.
plot_grid_km <- function(df, title, xlim, ylim, fill_limits) {
  ggplot(df, aes(X_km, Y_km, fill = Population)) +
    geom_tile() +
    coord_equal(xlim = xlim, ylim = ylim) +
    scale_fill_viridis_c(trans = "log10", limits = fill_limits) +
    labs(title = title, x = "km", y = "km") +
    theme_minimal()
}

for (key in names(countries)) {
  iso <- ISO_LOOKUP[[key]]
  real_csv <- file.path("data", iso, "Scales", "Raster550", "groups.csv")
  if (!file.exists(real_csv)) next

  real_groups  <- read_csv(real_csv, show_col_types = FALSE)
  synth_groups <- synth_results[[key]]$data

  real_c  <- to_km(real_groups, 550)     # real groups.csv is in metres
  synth_c <- synth_groups |> mutate(X_km = X, Y_km = Y)   # synthetic is already in km

  shared_xlim   <- range(c(real_c$X_km, synth_c$X_km))
  shared_ylim   <- range(c(real_c$Y_km, synth_c$Y_km))
  shared_limits <- range(c(real_c$Population, synth_c$Population))

  p <- plot_grid_km(real_c,  sprintf("%s Raster550 (real), n=%d, pop=%d", iso, nrow(real_c), sum(real_c$Population)),
                    shared_xlim, shared_ylim, shared_limits) +
       plot_grid_km(synth_c, sprintf("%s Synthetic, n=%d, pop=%d",       iso, nrow(synth_c), sum(synth_c$Population)),
                    shared_xlim, shared_ylim, shared_limits)
  show(p)
}

# --------------------------------------------------------------------------
# Diagnostic: the raster real_cell_populations() samples from, per country —
# 0 (real, uninhabited land) shown as a colour, not the NA (ocean) colour.
# --------------------------------------------------------------------------

for (key in names(countries)) {
  iso <- ISO_LOOKUP[[key]]
  r   <- real_population_raster(key, countries[[key]])
  pts <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  colnames(pts) <- c("x", "y", "Population")

  p <- ggplot(pts, aes(x, y, fill = Population)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(na.value = "grey85") +
    labs(title = sprintf("%s: raster sampled for Synth550 (pop=%d)", iso, round(sum(pts$Population, na.rm = TRUE))),
         x = NULL, y = NULL) +
    theme_minimal()
  show(p)
}
