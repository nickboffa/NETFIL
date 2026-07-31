library(tidyverse)
library(abc)
library(HDInterval)

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
message("Clearing stale population cache...")
system2("bash", args = c(file.path(model_dir, "clean_inputs.sh")), stdout = FALSE, stderr = FALSE)

# ── Observed targets ────────────────────────────────────────────────────────────
OBS_ANT_2016 <- 6.2   # ANT prevalence % (Lau et al. 2020, community survey age ≥8)
OBS_MF_2016  <- 1.59  # MF prevalence % (25.6% of ANT+ Mf-positive, back-calculated)
OBS_ICC_2016 <- 0.551603405400143 # linear from write_clustering_params.R for R550

# ── Semi-informative uniform priors (guessed from past plots)
T1_MIN <- 0.0;  T1_MAX <- 0.01
K_MIN  <- 0.0;  K_MAX  <- 0.3
W_MIN  <- 0.0;  W_MAX  <- 0.8

# ── ABC settings ───────────────────────────────────────────────────────────────
N_PARTICLES <- 500     # was 500
N_REPS      <- 2     # was 5
N_POSTERIOR <- 500   # was 10000, only affects output resolution, not fitted parameters
ABC_TOL     <- 0.05  # keep closest 5%

THETA2_VALS   <- c(1.0)          # just one to test
THETA2_LABELS <- c("Theta_1")
THETA2_SUFFIX <- c("1")


# ── Data setup ─────────────────────────────────────────────────────────────────
message("Copying ", COUNTRY, "/", SCALE, " scale data to data/...")
scale_files <- list.files(file.path(data_dir, COUNTRY, "Scales", SCALE), full.names = TRUE)
file.copy(scale_files, data_dir, overwrite = TRUE)
file.copy(file.path(data_dir, COUNTRY, "mda_params.csv"), data_dir, overwrite = TRUE)

# ── Helper: run model from model/ so DATADIR and OUTDIR resolve correctly ───────
run_model <- function(id) {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(model_dir)
  system2("./main", args = id, stdout = FALSE, stderr = FALSE)
}

# ── Helper: parse output CSV for year-start summary statistics ─────────────────
# output_epidemics writes quarterly; we use year == 2016, day == 0 rows.
# Returns one row per simulation replicate with columns ANT, MF, ICC.
parse_output <- function(path) {
  tryCatch({
    # `name` is a string column; everything else is numeric
    dat <- read_csv(path, col_types = cols(name = "c", .default = "d"), show_col_types = FALSE)
    r16 <- dplyr::filter(dat, year == 2016, day == 0)
    if (nrow(r16) == 0) return(tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_))

    # ICC from per-group MF prevalences: ICC = sigma2 / (sigma2 + pi^2/3)
    # where sigma2 = var(logit(p_j)) across groups
    calc_icc <- function(row) {
      mf_vals  <- unlist(dplyr::select(row, matches("^mf_\\d")))
      pop_vals <- unlist(dplyr::select(row, matches("^pop_\\d")))
      p_j <- mf_vals / pop_vals
      p_j <- p_j[is.finite(p_j) & p_j > 0 & p_j < 1]
      if (length(p_j) < 2) return(NA_real_)
      sigma2 <- var(log(p_j / (1 - p_j)))
      sigma2 / (sigma2 + pi^2 / 3)
    }

    tibble(
      ANT = if_else(r16$pop_total > 0, r16$ant_total / r16$pop_total * 100, NA_real_),
      MF  = if_else(r16$pop_total > 0, r16$mf_total  / r16$pop_total * 100, NA_real_),
      ICC = purrr::map_dbl(seq_len(nrow(r16)), \(i) calc_icc(r16[i, ]))
    )
  }, error = function(e) {
    warning("Failed to parse ", path, ": ", conditionMessage(e))
    tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_)
  })
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
  colnames(post) <- c("T1", "W", "k")

  theta2_val <- THETA2_VALS[THETA2_LABELS == label]
  write_csv(tibble(
    Theta_1   = median(post$T1),
    Theta_2   = theta2_val,
    Agg       = median(post$k),
    WorktoNot = median(post$W)
  ), file.path(out_dir, "tran_params.csv"))

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
  adj_df <- as.data.frame(abc_fit$adj.values) |> setNames(c("T1", "W", "k"))
  save_png(
    bind_rows(
      prior  |> mutate(type = "Prior"),
      adj_df |> mutate(type = "Posterior")
    ) |>
      pivot_longer(c(T1, W, k), names_to = "param", values_to = "value") |>
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

  message(sprintf("\n── %s  %s (theta2 = %s) ──────────────────────────",
                  SCALE, label, theta2))

  csv_name <- sprintf("fit_%s.csv", suffix)
  csv_path <- file.path(out_theta, csv_name)
  if (file.exists(csv_path)) file.remove(csv_path)  # start fresh each run

  results_list <- vector("list", N_PARTICLES)

  for (p in seq_len(N_PARTICLES)) {

    t1 <- runif(1, T1_MIN, T1_MAX)
    k  <- runif(1, K_MIN,  K_MAX)
    w  <- runif(1, W_MIN,  W_MAX)

    write_csv(
      tibble(Theta_1 = t1, Theta_2 = theta2, Agg = k, WorktoNot = w),
      file.path(data_dir, "tran_params.csv")
    )

    id       <- sprintf("%s_p%04d", suffix, p)
    out_file <- file.path(output_dir, id)

    message(sprintf("  Running particle %d / %d", p, N_PARTICLES))
    run_model(id)

    if (file.exists(out_file)) {
      stats <- parse_output(out_file)
      file.remove(out_file)
      results_list[[p]] <- mutate(stats, Sim = p, T1 = t1, W = w, k = k)
    } else {
      message(sprintf("  Warning: no output for p=%d", p))
    }
    
    message(sprintf("  %d / %d particles", p, N_PARTICLES))
    if (p %% 5 == 0) {
      write_csv(bind_rows(results_list), csv_path)
    }
  }

  particles <- bind_rows(results_list)

  # ── ABC-GLM ──────────────────────────────────────────────────────────────────
  complete <- particles |>
    drop_na(ANT, MF, ICC) |>
    group_by(Sim, T1, W, k) |>
    summarise(
      ANT = mean(ANT),
      MF  = mean(MF),
      ICC = mean(ICC),
      .groups = "drop"
    )
  
  message(sprintf("  %d / %d complete rows for ABC", nrow(complete), nrow(particles)))
  
  # L1 rejection: keep closest ABC_TOL fraction by L1 distance
  obs_vec <- c(OBS_ANT_2016, OBS_MF_2016, OBS_ICC_2016)
  ss_mat  <- as.matrix(dplyr::select(complete, ANT, MF, ICC))
  l1_dist <- rowSums(abs(sweep(ss_mat, 2, obs_vec)))
  accepted <- complete[l1_dist <= quantile(l1_dist, ABC_TOL), ]
  message(sprintf("  %d particles accepted after L1 rejection (tol = %.0f%%)",
                  nrow(accepted), ABC_TOL * 100))

  # GLM post-sampling adjustment on accepted particles
  fit <- abc(
    target  = obs_vec,
    param   = dplyr::select(accepted, T1, W, k),
    sumstat = dplyr::select(accepted, ANT, MF, ICC),
    tol     = 1.0,
    method  = "loclinear",
    transf  = c("logit", "logit", "logit"),
    logit.bounds = matrix(c(T1_MIN, T1_MAX,
                            W_MIN,  W_MAX,
                            K_MIN,  K_MAX),
                          nrow = 3, byrow = TRUE)
  )

  prior <- tibble(
    T1 = runif(N_POSTERIOR, T1_MIN, T1_MAX),
    W  = runif(N_POSTERIOR, W_MIN,  W_MAX),
    k  = runif(N_POSTERIOR, K_MIN,  K_MAX)
  )
  write_abc_glm_files(fit, out_theta, label, prior,
                      complete, accepted, l1_dist, obs_vec)
  message(sprintf("  Written to %s", out_theta))
}

message("\nDone. Results in ", fit_base)

variation_dat |> 
  filter(is.na(value))


