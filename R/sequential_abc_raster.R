library(tidyverse)
library(abc)        # still used for the final GLM regression-adjustment step
library(EasyABC)    # ABC_sequential() — sequential ABC / ABC-SMC implementation
library(HDInterval)
library(jsonlite)

# ── Country / scale configuration ─────────────────────────────────────────────
COUNTRY <- "ASM"
SCALE   <- "Raster550"

# ── Paths ──────────────────────────────────────────────────────────────────────
script_dir <- here::here()
model_dir  <- file.path(script_dir, "model")
data_dir   <- file.path(script_dir, "data")
output_dir <- file.path(script_dir, "output")   # = model/../output = ../output from binary
fit_base   <- file.path(data_dir, COUNTRY, "Fitted", SCALE)
binary     <- file.path(model_dir, "main")

dir.create(output_dir, showWarnings = FALSE)
dir.create(fit_base,   showWarnings = FALSE, recursive = TRUE)

# ── Clear stale population state ───────────────────────────────────────────────
# ABC.init from a previous scale would cause a crash (wrong group_blocks loaded).
# clean_inputs.sh uses paths relative to model/ (../$config/...) — it must be
# run with model/ as the working directory, or those paths silently resolve
# to the wrong place and nothing gets deleted (this WAS happening: relative
# to script_dir the "../$config/" landed one level above the repo entirely).
message("Clearing stale population cache...")
old_wd <- setwd(model_dir)
system2("bash", args = "clean_inputs.sh", stdout = FALSE, stderr = FALSE)
setwd(old_wd)

# ── Observed targets ────────────────────────────────────────────────────────────
OBS_ANT_2016 <- 6.2   # ANT prevalence % (Lau et al. 2020, community survey age ≥8)
OBS_MF_2016  <- 1.59  # MF prevalence % (25.6% of ANT+ Mf-positive, back-calculated)
OBS_ICC_2016 <- 0.49 # 0.5134347 # exponential y = a * exp (b * x)

#OBS_ICC_2016 <- 0.551603405400143 # village-level ICC estimated from 2016 survey's
                                   # household- and region-level ICC (see notes) —
                                   # independent of the 2010 seeding ICC below, which
                                   # is now fitted rather than fixed by distance decay

# ICC is now load-bearing: it's the only thing constraining the fitted 2010
# seeding ICC (see ICC10_MIN/MAX below), so it must be fit. Kept as a toggle
# only for quick debugging — leave TRUE for real fits.
FIT_ICC <- TRUE

# ── Semi-informative uniform priors (guessed from past plots)
T1_MIN <- 0.0;  T1_MAX <- 0.015 # 0.01 originally
K_MIN  <- 0.01;  K_MAX  <- 0.75 # 0.3 originally
W_MIN  <- 0.0;  W_MAX  <- 0.8

# T1_MIN <- 0.0;    T1_MAX <- 0.03
# K_MIN  <- 0.0;  K_MAX  <- 0.6
# W_MIN  <- 0.0;    W_MAX  <- 0.2

# k -> 0 makes bite_gamma(0, Inf) hang forever during agent construction, so
# the prior actually sampled from below floors K_MIN at this value.
K_MIN_EFF <- max(K_MIN, 0.01)

# 2010 seeding ICC — previously fixed via the distance-decay line in
# write_clustering_params.R; now a free parameter fit against OBS_ICC_2016.
# ICC10_MIN <- 0.1;  ICC10_MAX <- 0.6
ICC10_MIN <- 0.4; ICC10_MAX <- 0.9


# ── ABC settings ───────────────────────────────────────────────────────────────
N_PARTICLES <- 50      # size of the final retained population — see cost note below
N_REPS      <- 2     # was 5
N_POSTERIOR <- 10000   # was 10000, only affects output resolution, not fitted parameters
ABC_TOL     <- 0.05  # only used cosmetically now, for the L1-distance plot's reference line
N_WORKERS   <- 6     # parallel workers (= cores to use; leave ≥1 free for OS/VSCode)

# ── Sequential ABC (ABC-SMC) settings ───────────────────────────────────────────
# EasyABC::ABC_sequential() replaces abc_raster.R's one-shot "draw N_PARTICLES,
# keep the closest ABC_TOL fraction" step. It runs Lenormand et al. (2012)'s
# adaptive population Monte Carlo algorithm: an initial population is drawn
# from the prior, then repeatedly resampled + perturbed + re-filtered against
# a shrinking tolerance until the per-step acceptance rate drops below
# SMC_P_ACC_MIN. This is EasyABC's own recommended default method — of its
# five algorithms ("Beaumont", "Drovandi", "Delmoral", "Lenormand",
# "Emulation") it's the one built specifically for uniform priors (which is
# all we have here) and needs no hand-tuned tolerance schedule.
#
# `nb_simul` is EasyABC's total per-step simulation budget, of which
# `alpha` (kept) survive each step as the retained population — so to end up
# with a final population of N_PARTICLES, nb_simul = N_PARTICLES / SMC_ALPHA.
# Total simulation cost is nb_simul for the first step, plus
# nb_simul * (1 - SMC_ALPHA) for every step after, until p_acc < SMC_P_ACC_MIN.
#
# COST NOTE: sequential ABC only costs *less* than abc_raster.R's one-shot
# N_PARTICLES=600 draw if it's not also asked for a much larger/tighter final
# population — abc_raster.R actually only ever kept ~ABC_TOL*600 = 30
# particles as "accepted". Asking EasyABC for a final population of 600 (an
# earlier version of this script did) meant step 1 alone already cost 1200,
# before any adaptive refinement steps — several times abc_raster.R's total.
# With N_PARTICLES=50 and SMC_ALPHA=0.5: step 1 costs 100, and every
# subsequent step costs another 50, stopping once a step's acceptance rate
# drops below SMC_P_ACC_MIN. Raising SMC_P_ACC_MIN (loosening the stopping
# rule) trades some final-tolerance sharpness for fewer of those extra steps
# — EasyABC's own docs frame it exactly this way. With these settings, total
# cost should land at a few hundred simulations — comparably fast to, or
# faster than, abc_raster.R's 600 — though (unlike a fixed-generation-count
# design) the exact total isn't knowable in advance: it depends on how
# quickly the acceptance rate falls, which is problem-dependent.
SMC_ALPHA     <- 0.5   # EasyABC default: fraction of each step's draws kept
SMC_P_ACC_MIN <- 0.15  # raised from EasyABC's 0.05 default — stop sooner, trade some precision for speed

THETA2_VALS   <- c(1.0)          # just one to test
THETA2_LABELS <- c("Theta_1")
THETA2_SUFFIX <- c("1")


# ── Data setup ─────────────────────────────────────────────────────────────────
message("Copying ", COUNTRY, "/", SCALE, " scale data to data/...")
scale_files <- list.files(file.path(data_dir, COUNTRY, "Scales", SCALE), full.names = TRUE)
file.copy(scale_files, data_dir, overwrite = TRUE)
file.copy(file.path(data_dir, COUNTRY, "country.json"), data_dir, overwrite = TRUE)
writeLines(
  c(paste0("country=", COUNTRY), paste0("scale=", SCALE), paste0("theta2=", THETA2_SUFFIX[1])),
  file.path(data_dir, "current_state.txt")
)

# Target 2010 seeded ANT prevalence — needed to solve beta_0 for each
# particle's sampled ICC10 (see icc_to_sigma_beta below).
country_cfg <- fromJSON(file.path(data_dir, COUNTRY, "country.json"))
ANT_0 <- country_cfg$seeding$init_ant_prev

# ── Helper: convert a 2010 seeding ICC into (sigma_group, beta_0) ──────────────
# Same derivation as write_clustering_params.R's calc_sigma_and_beta(), but
# taking ICC directly instead of computing it from inter-group distance —
# ICC10 is now a free ABC parameter, not fixed by the distance-decay line.
# Defined at top level (not just inside build_model() below) because the
# warm-up run outside the parallel workers needs it too.
icc_to_sigma_beta <- function(icc, ant_0) {
  sigma_group <- sqrt(icc / (1 - icc) * pi^2 / 3)

  expected_prev <- function(mu) {
    integrate(
      function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma_group),
      lower = -Inf, upper = Inf
    )$value
  }

  # expected_prev(mu) is monotonic increasing in mu, from ~0 (mu -> -Inf) to
  # 0.5 (mu = 0), so any sufficiently negative lower bound brackets a root.
  # A fixed -10 worked for the old ICC10_MAX (0.6, sigma_group ~2.2) but not
  # the wider prior now in use: at ICC10 = 0.9, sigma_group ~5.4 and
  # expected_prev(-10) is already above target ant_0, so uniroot() throws
  # "f() values at end points not of opposite sign" for high-ICC10 particles
  # (intermittent, since ICC10 is drawn per particle). Scale the lower bound
  # with sigma_group so it stays valid regardless of how wide the ICC10 prior
  # gets.
  beta_0 <- uniroot(
    f        = function(mu) expected_prev(mu) - ant_0,
    interval = c(-10 - 10 * sigma_group, 0)
  )$root

  list(sigma_group = sigma_group, beta_0 = beta_0)
}

# ── Helper: build the model() function ABC_sequential() calls per particle ─────
# EasyABC's cluster mode (n_cluster > 1) dispatches `model` to PSOCK workers
# via parallel::parLapply *without* exporting the calling session's globals —
# each worker is a bare fresh R process (only base/methods/datasets/utils/
# stats attached, no tidyverse). So everything `model` needs — model_dir,
# output_dir, theta2, suffix, N_REPS, ANT_0, fit_stats, and the run_model/
# parse_output helpers themselves — has to live in model's own closure
# (captured here as build_model()'s local variables) rather than being read
# off the top-level script environment, and every tidyverse call inside has
# to be fully namespaced (readr::, dplyr::, purrr::) since nothing is
# library()'d on the worker.
#
# use_seed = TRUE (required by EasyABC for n_cluster > 1 anyway) hands model()
# a per-call integer seed as x[1] — used here not for RNG seeding (the model's
# randomness lives in the C++ binary, not R) but as a collision-free id for
# run_model()'s output file, unique across every particle in the whole
# sequential run.
build_model <- function(model_dir, output_dir, theta2, suffix, N_REPS, ANT_0, fit_stats) {
  # These arguments are lazy promises pointing at symbols in the *master*
  # session's global environment. If left unforced, that lookup is still
  # unresolved when this closure is serialized to a worker — and a promise's
  # environment doesn't survive serialization the way a plain binding does
  # (globalenv is special-cased and comes back empty on the worker). Forcing
  # every argument here, before it's captured by the nested functions below,
  # bakes in its actual value instead of a dangling lookup.
  force(model_dir); force(output_dir); force(theta2); force(suffix)
  force(N_REPS); force(ANT_0); force(fit_stats)

  # Same derivation as the top-level icc_to_sigma_beta() (used by the warm-up
  # run and write_abc_glm_files()) — duplicated here, rather than called by
  # name, so it's part of *this* closure's own environment and travels with
  # it to the workers. A reference to the top-level version would suffer the
  # same unresolved-lookup problem described above, since function lookups
  # are resolved lexically through this frame's parent (globalenv) too.
  icc_to_sigma_beta <- function(icc, ant_0) {
    sigma_group <- sqrt(icc / (1 - icc) * pi^2 / 3)

    expected_prev <- function(mu) {
      integrate(
        function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma_group),
        lower = -Inf, upper = Inf
      )$value
    }

    beta_0 <- uniroot(
      f        = function(mu) expected_prev(mu) - ant_0,
      interval = c(-10 - 10 * sigma_group, 0)
    )$root

    list(sigma_group = sigma_group, beta_0 = beta_0)
  }

  run_model <- function(id, t1, t2, k, w, sigma_group, beta_0) {
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(model_dir)
    # output_epidemics() appends (ios::app) to this file — if a leftover file from
    # a crashed prior run with the same id is still sitting in output/, the new
    # run's rows get silently mixed in with stale ones. Remove it up front so a
    # crash never leaves contamination for a later run to inherit.
    unlink(file.path(output_dir, id))
    system2("./main", args = c(id, N_REPS, t1, t2, k, w, sigma_group, beta_0), stdout = FALSE, stderr = FALSE)
  }

  # Parses output CSV for year-start summary statistics. output_epidemics
  # writes quarterly; we use year == 2016, day == 0 rows. Returns one row per
  # simulation replicate with columns ANT, MF, ICC.
  parse_output <- function(path) {
    tryCatch({
      # `name` is a string column; everything else is numeric
      dat <- readr::read_csv(path, col_types = readr::cols(name = "c", .default = "d"), show_col_types = FALSE)
      r16 <- dplyr::filter(dat, year == 2016, day == 0)
      if (nrow(r16) == 0) return(tibble::tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_))

      # ICC from per-group MF prevalences: ICC = sigma2 / (sigma2 + pi^2/3)
      # where sigma2 = var(logit(p_j)) across groups
      calc_icc <- function(row) {
        mf_vals  <- unlist(dplyr::select(row, dplyr::matches("^mf_\\d")))
        pop_vals <- unlist(dplyr::select(row, dplyr::matches("^pop_\\d")))
        p_j <- mf_vals / pop_vals
        p_j <- p_j[is.finite(p_j) & p_j > 0 & p_j < 1]
        if (length(p_j) < 2) return(NA_real_)
        sigma2 <- var(log(p_j / (1 - p_j)))
        sigma2 / (sigma2 + pi^2 / 3)
      }

      tibble::tibble(
        ANT = dplyr::if_else(r16$pop_total > 0, r16$ant_total / r16$pop_total * 100, NA_real_),
        MF  = dplyr::if_else(r16$pop_total > 0, r16$mf_total  / r16$pop_total * 100, NA_real_),
        ICC = purrr::map_dbl(seq_len(nrow(r16)), \(i) calc_icc(r16[i, ]))
      )
    }, error = function(e) {
      warning("Failed to parse ", path, ": ", conditionMessage(e))
      tibble::tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_)
    })
  }

  function(x) {
    seed  <- as.integer(round(x[1]))
    t1    <- x[2]; w <- x[3]; k <- x[4]; icc10 <- x[5]
    cb <- icc_to_sigma_beta(icc10, ANT_0)

    id <- sprintf("%s_s%d", suffix, seed)
    run_model(id, t1, theta2, k, w, cb$sigma_group, cb$beta_0)

    out_file <- file.path(output_dir, id)
    if (!file.exists(out_file)) {
      # Model failed to produce output (rare). ABC_sequential needs a numeric
      # vector back from every call, so return something far from any
      # plausible target rather than NA — the particle is then naturally
      # rejected by the distance filter instead of crashing the sampler.
      return(rep(1e6, length(fit_stats)))
    }
    stats <- parse_output(out_file)
    file.remove(out_file)

    icc_vals <- stats$ICC
    icc_vals[is.na(icc_vals)] <- 0  # occurs when MF prevalence is really low
    result <- c(ANT = mean(stats$ANT), MF = mean(stats$MF), ICC = mean(icc_vals))
    as.numeric(result[fit_stats])
  }
}

# ── Helper: write ABC-GLM output files ─────────────────────────────────────────
# Writes the minimum needed to reload results (abc_fit.rds + posterior draws +
# TranParams.csv) plus key diagnostic PNGs.
write_abc_glm_files <- function(abc_fit, out_dir, label, prior,
                                 particles, accepted, l1_dist, obs_vec) {
  # ── Reload essentials ────────────────────────────────────────────────────────
  write_rds(abc_fit, file.path(out_dir, "abc_fit.rds"))

  adj <- as.data.frame(abc_fit$adj.values)
  wts <- abc_fit$weights / sum(abc_fit$weights)
  idx  <- sample(nrow(adj), N_POSTERIOR, replace = TRUE, prob = wts)
  post <- adj[idx, ]
  colnames(post) <- c("T1", "W", "k", "ICC10")

  theta2_val <- THETA2_VALS[THETA2_LABELS == label]
  write_csv(tibble(
    Theta_1   = median(post$T1),
    Theta_2   = theta2_val,
    Agg       = median(post$k),
    WorktoNot = median(post$W),
    ICC10 = median(post$ICC10) # isn't actually needed because of clustering_params.csv
  ), file.path(out_dir, "tran_params.csv"))

  # Fitted clustering params — convert the posterior median ICC10 back into
  # (sigma_group, beta_0) so a production (non-ABC) run can use it directly.
  # The GLM local-linear adjustment (method="loclinear") doesn't respect
  # parameter bounds — with few accepted particles it can push the adjusted
  # median outside (0, 1), where sigma_group = sqrt(icc/(1-icc) * pi^2/3) and
  # the beta_0 root-find blow up (observed directly while testing this).
  # Clamp defensively rather than crash.
  icc10_median <- median(post$ICC10)
  if (icc10_median <= 0 || icc10_median >= 1) {
    warning(sprintf(
      "Posterior median ICC10 (%.4f) is outside (0, 1) — GLM adjustment likely extrapolated past valid bounds (check accepted-particle count). Clamping to fit.",
      icc10_median
    ))
    icc10_median <- min(max(icc10_median, 0.001), 0.999)
  }
  fitted_cluster <- icc_to_sigma_beta(icc10_median, ANT_0)
  write_csv(tibble(
    sigma_group = fitted_cluster$sigma_group,
    beta_0      = fitted_cluster$beta_0
  ), file.path(out_dir, "clustering_params.csv"))

  # ── Diagnostic PNGs ──────────────────────────────────────────────────────────
  save_png <- function(p, name, w = 10, h = 6) {
    ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 150)
  }

  # Inter-run vs overall variation
  variation_dat <- particles |>
    pivot_longer(c(ANT, MF, ICC), names_to = "stat", values_to = "value")

  noise_labels <- variation_dat |>
    group_by(stat) |>
    mutate(total_sd = sd(value, na.rm = TRUE)) |>
    group_by(stat, Sim) |>
    summarise(within_var = var(value, na.rm = TRUE), total_sd = first(total_sd),
              .groups = "drop") |>
    group_by(stat) |>
    summarise(
      noise_frac = sqrt(mean(within_var, na.rm = TRUE)) / first(total_sd),
      .groups = "drop"
    ) |>
    mutate(label = sprintf("%s (noise: %.0f%%)", stat, noise_frac * 100)) |>
    dplyr::select(stat, label) |>
    deframe()

  var_plot <- ggplot(variation_dat, aes(x = Sim, y = value, color = factor(Sim))) +
    geom_point(alpha = 0.5, size = 0.8) +
    facet_wrap(~stat, scales = "free_y", ncol = 1,
               labeller = as_labeller(noise_labels)) +
    labs(x = "Particle", y = NULL) +
    guides(color = "none") +
    theme_minimal()

  save_png(
    var_plot,
    "variation.png", h = 10
  )

  # L1 distance histogram with acceptance cutoff
  save_png(
    ggplot(tibble(l1 = l1_dist), aes(x = l1)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "white") +
      geom_vline(xintercept = quantile(l1_dist, ABC_TOL),
                 color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = quantile(l1_dist, ABC_TOL), y = Inf,
               label = sprintf("%.0f%% cutoff", ABC_TOL * 100),
               hjust = -0.1, vjust = 2, color = "red") +
      labs(title = "L1 distance distribution",
           x = "L1 distance to observed", y = "Count") +
      theme_minimal(),
    "l1_distances.png"
  )

  # All particles vs accepted vs observed (two facets: MF and ICC on y-axis)
  obs_labels <- c(MF = "MF prevalence 2016 (%)", ICC = "ICC")
  particle_long <- bind_rows(
    particles |> mutate(y = MF,  facet = "MF"),
    particles |> mutate(y = ICC, facet = "ICC")
  )
  accepted_long <- bind_rows(
    accepted  |> mutate(y = MF,  facet = "MF"),
    accepted  |> mutate(y = ICC, facet = "ICC")
  )
  obs_long <- tibble(
    x     = c(obs_vec[1], obs_vec[1]),
    y     = c(obs_vec[2], obs_vec[3]),
    facet = c("MF", "ICC")
  )
  save_png(
    ggplot() +
      geom_point(data = particle_long, aes(ANT, y),
                 color = "grey70", alpha = 0.4, size = 1) +
      geom_point(data = accepted_long, aes(ANT, y),
                 color = "steelblue", alpha = 0.6, size = 1.5) +
      geom_point(data = obs_long, aes(x, y),
                 color = "red", shape = 8, size = 4, stroke = 1.5) +
      facet_wrap(~facet, scales = "free_y",
                 labeller = as_labeller(obs_labels)) +
      labs(title = "All particles (grey) vs accepted (blue) vs observed (red)",
           x = "Antigen prevalence 2016 (%)") +
      theme_minimal() +
      theme(axis.title.y = element_blank()),
    "particles.png", w = 14
  )

  # Prior vs posterior (regression-adjusted)
  adj_df <- as.data.frame(abc_fit$adj.values) |> setNames(c("T1", "W", "k", "ICC10"))
  save_png(
    bind_rows(
      prior  |> mutate(type = "Prior"),
      adj_df |> mutate(type = "Posterior")
    ) |>
      pivot_longer(c(T1, W, k, ICC10), names_to = "param", values_to = "value") |>
      ggplot(aes(x = value, color = type, fill = type)) +
      geom_density(alpha = 0.2) +
      facet_wrap(~param, scales = "free") +
      scale_color_manual(values = c("Prior" = "grey50", "Posterior" = "firebrick")) +
      scale_fill_manual(values  = c("Prior" = "grey50", "Posterior" = "firebrick")) +
      labs(title = "Prior vs posterior (regression-adjusted)",
           x = "Value", y = "Density", color = NULL, fill = NULL) +
      theme_minimal(),
    "posterior.png", w = 12, h = 5
  )
}

# ── Main: loop over theta2 values ─────────────────────────────────────────────
for (i in seq_along(THETA2_VALS)) {

  theta2 <- THETA2_VALS[i]
  label  <- THETA2_LABELS[i]
  suffix <- THETA2_SUFFIX[i]
  out_theta <- file.path(fit_base, label)
  dir.create(out_theta, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n── %s  %s (theta2 = %s) — sequential ABC (EasyABC::ABC_sequential, Lenormand) ──",
                  SCALE, label, theta2))

  csv_name <- sprintf("fit_%s.csv", suffix)
  csv_path <- file.path(out_theta, csv_name)
  if (file.exists(csv_path)) file.remove(csv_path)  # start fresh each run

  # Warm-up: pre-build the .init population cache so parallel workers don't race to
  # create it once ABC_sequential() starts dispatching to its cluster. K_MIN can be 0,
  # which causes bite_gamma(0, Inf) to hang on agent births — clamp to K_MIN_EFF.
  # Clustering params don't affect population construction, so any valid value works
  # here — use the midpoint of the ICC10 prior. Population is independent of
  # (T1, k, W, ICC10), so this only needs to happen once per theta2.
  message("  Warm-up run to build population cache...")
  warmup_cluster <- icc_to_sigma_beta((ICC10_MIN + ICC10_MAX) / 2, ANT_0)
  t_warmup_start <- Sys.time()
  { # same body as build_model()'s run_model(), run once from the top level
    old_wd <- setwd(model_dir)
    unlink(file.path(output_dir, "warmup"))
    system2("./main", args = c("warmup", N_REPS, max(T1_MIN, 1e-4), theta2, K_MIN_EFF,
                                max(W_MIN, 1e-4), warmup_cluster$sigma_group, warmup_cluster$beta_0),
            stdout = FALSE, stderr = FALSE)
    setwd(old_wd)
  }
  t_warmup_end <- Sys.time()
  suppressWarnings(file.remove(file.path(output_dir, "warmup")))
  warmup_sec <- as.numeric(difftime(t_warmup_end, t_warmup_start, units = "secs"))

  fit_stats <- if (FIT_ICC) c("ANT", "MF", "ICC") else c("ANT", "MF")
  obs_fit   <- if (FIT_ICC) c(OBS_ANT_2016, OBS_MF_2016, OBS_ICC_2016) else
               c(OBS_ANT_2016, OBS_MF_2016)
  obs_vec   <- c(OBS_ANT_2016, OBS_MF_2016, OBS_ICC_2016)  # always 3-element for visualisation

  # ── Sequential ABC ────────────────────────────────────────────────────────────
  prior <- list(
    c("unif", T1_MIN, T1_MAX),
    c("unif", W_MIN,  W_MAX),
    c("unif", K_MIN_EFF, K_MAX),
    c("unif", ICC10_MIN, ICC10_MAX)
  )
  model <- build_model(model_dir, output_dir, theta2, suffix, N_REPS, ANT_0, fit_stats)

  nb_simul_arg <- ceiling(N_PARTICLES / SMC_ALPHA)
  message(sprintf("  Running ABC_sequential (target %d particles/step, alpha = %.2f, p_acc_min = %.2f, %d workers)...",
                  N_PARTICLES, SMC_ALPHA, SMC_P_ACC_MIN, N_WORKERS))
  t_smc_start <- Sys.time()
  smc <- ABC_sequential(
    method              = "Lenormand",
    model               = model,
    prior               = prior,
    nb_simul            = nb_simul_arg,
    summary_stat_target = obs_fit,
    use_seed            = TRUE,
    n_cluster           = N_WORKERS,
    alpha               = SMC_ALPHA,
    p_acc_min           = SMC_P_ACC_MIN,
    verbose = TRUE
  )
  t_smc_end <- Sys.time()
  smc_sec <- as.numeric(difftime(t_smc_end, t_smc_start, units = "secs"))
  message(sprintf("  Wall time: %.0fs | %d total simulations | final tolerance (epsilon) = %.4f",
                  smc_sec, smc$nsim, smc$epsilon))

  write_csv(
    tibble(
      scale             = SCALE,
      theta2            = theta2,
      method             = "Lenormand (EasyABC::ABC_sequential)",
      n_particles_final = nrow(smc$param),
      alpha             = SMC_ALPHA,
      p_acc_min         = SMC_P_ACC_MIN,
      n_reps            = N_REPS,
      n_workers         = N_WORKERS,
      n_sim_total       = smc$nsim,
      epsilon           = smc$epsilon,
      warmup_sec        = round(warmup_sec, 1),
      wall_sec          = round(smc_sec, 1),
      package_computime_sec = round(smc$computime, 1),
      started_at        = format(t_smc_start, "%Y-%m-%d %H:%M:%S"),
      finished_at       = format(t_smc_end,   "%Y-%m-%d %H:%M:%S")
    ),
    file.path(out_theta, "timing.csv")
  )

  # ── Build the accepted-particle tibble from EasyABC's final population ────────
  # ABC_sequential()'s final population *is* its accepted set (by construction —
  # every returned particle already passed the algorithm's own final-step
  # tolerance), unlike abc_raster.R's separate "draw then reject" split. So the
  # "all particles" and "accepted" views below are the same population; passed
  # to write_abc_glm_files() twice for its grey-vs-blue plot to keep that
  # function unchanged.
  colnames(smc$param) <- c("T1", "W", "k", "ICC10")
  colnames(smc$stats)  <- fit_stats

  complete <- as_tibble(smc$param) |>
    bind_cols(as_tibble(smc$stats)) |>
    mutate(Sim = row_number(), weight = smc$weights)
  if (!FIT_ICC) complete$ICC <- NA_real_  # write_abc_glm_files always plots ICC

  ss_mat  <- as.matrix(dplyr::select(complete, all_of(fit_stats)))
  l1_dist <- rowSums(abs(sweep(ss_mat, 2, obs_fit)))

  message(sprintf("  %d particles in final population", nrow(complete)))

  # ── GLM post-sampling adjustment ───────────────────────────────────────────────
  # Resample the final weighted population down to an unweighted particle set —
  # particles that carried more SMC weight appear more often — so it drops
  # straight into the same local-linear regression-adjustment step
  # abc_raster.R uses.
  resample_idx <- sample(nrow(complete), nrow(complete), replace = TRUE, prob = complete$weight)
  accepted     <- complete[resample_idx, ]

  fit <- abc(
    target  = obs_fit,
    param   = dplyr::select(accepted, T1, W, k, ICC10),
    sumstat = dplyr::select(accepted, all_of(fit_stats)),
    tol     = 1.0,
    method  = "loclinear",
    transf  = "none"
  )

  prior_draws <- tibble(
    T1    = runif(N_POSTERIOR, T1_MIN, T1_MAX),
    W     = runif(N_POSTERIOR, W_MIN,  W_MAX),
    k     = runif(N_POSTERIOR, K_MIN,  K_MAX),
    ICC10 = runif(N_POSTERIOR, ICC10_MIN, ICC10_MAX)
  )
  write_csv(complete, csv_path)
  write_abc_glm_files(fit, out_theta, label, prior_draws,
                      complete, accepted, l1_dist, obs_vec)

  message(sprintf("  Written to %s", out_theta))
}

message("\nDone. Results in ", fit_base)
