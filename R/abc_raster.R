library(tidyverse)
library(abc)
library(HDInterval)

# ── Scale configuration ────────────────────────────────────────────────────────
# Change SCALE to rerun for a different spatial scale.
#   "Raster660"  → fit_r*.tsv   (euc_dists.csv, no road_dist)
#   "Many"       → fit_s*.tsv   (uses existing euc_dist.csv + road_dist.csv)
#   "Village"    → fit_v*.tsv   (uses existing euc_dist.csv + road_dist.csv)
SCALE <- "Raster660"

TSV_PREFIX <- switch(SCALE,
  Raster660 = "r",
  Many      = "s",
  Village   = "v",
  stop("Unknown scale: ", SCALE)
)

# ── Paths ──────────────────────────────────────────────────────────────────────
script_dir <- here::here()
model_dir  <- file.path(script_dir, "model")
data_dir   <- file.path(script_dir, "data")
output_dir <- file.path(script_dir, "output")   # = model/../output = ../output from binary
fit_base   <- file.path(data_dir, "Fitted", SCALE)
binary     <- file.path(model_dir, "main")

dir.create(output_dir, showWarnings = FALSE)
dir.create(fit_base,   showWarnings = FALSE, recursive = TRUE)

# ── Clear stale population state ───────────────────────────────────────────────
# ABC.init from a previous scale would cause a crash (wrong group_blocks loaded).
message("Clearing stale population cache...")
system2("bash", args = c(file.path(model_dir, "clean_inputs.sh")), stdout = FALSE, stderr = FALSE)

# ── Observed targets ────────────────────────────────────────────────────────────
OBS_ANT_2016 <- 6.2   # Ag prevalence % (Lau et al. 2020, community survey age ≥8)
OBS_MF_2016  <- 1.59  # MF prevalence % (25.6% of Ag+ Mf-positive, back-calculated)

# ── Priors (semi-informative uniform, matching paper) ──────────────────────────
T1_MIN <- 0.0;  T1_MAX <- 0.01
K_MIN  <- 0.0;  K_MAX  <- 0.3
W_MIN  <- 0.0;  W_MAX  <- 0.8

# ── ABC settings ───────────────────────────────────────────────────────────────
N_PARTICLES <- 500     # was 500
N_REPS      <- 5     # was 5
N_POSTERIOR <- 2500   # was 10000 (needs to be <= accepted particles)
ABC_TOL     <- 0.05  # keep closest 5%

THETA2_VALS   <- c(1.0)          # just one to test
THETA2_LABELS <- c("Theta_1")
THETA2_SUFFIX <- c("1")


# ── Data setup ─────────────────────────────────────────────────────────────────
message("Copying ", SCALE, " scale data to data/...")
scale_dir <- file.path(data_dir, "Scales", SCALE)

file.copy(file.path(scale_dir, "groups.csv"), file.path(data_dir, "groups.csv"), overwrite = TRUE)

if (SCALE == "Raster660") {
  file.copy(file.path(scale_dir, "euc_dist.csv"), file.path(data_dir, "euc_dist.csv"), overwrite = TRUE)
  file.copy(file.path(scale_dir, "euc_dist.csv"), file.path(data_dir, "road_dist.csv"), overwrite = TRUE)
} else {
  file.copy(file.path(scale_dir, "euc_dist.csv"),  file.path(data_dir, "euc_dist.csv"),  overwrite = TRUE)
  file.copy(file.path(scale_dir, "road_dist.csv"), file.path(data_dir, "road_dist.csv"), overwrite = TRUE)
}

# ── Helper: run model from model/ so DATADIR and OUTDIR resolve correctly ───────
run_model <- function(id) {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(model_dir)
  system2("./main", args = id, stdout = FALSE, stderr = FALSE)
}

# ── Helper: parse output CSV for year-start summary statistics ─────────────────
# output_epidemics writes quarterly; we use Day == 0 rows.
# Column positions (1-indexed after read_csv): 2=Year 3=Day 21=Pop 22=inf 23=ant
parse_output <- function(path) {
  tryCatch({
    dat <- read_csv(path, col_types = cols(.default = "d"), show_col_types = FALSE)

    r14 <- dplyr::filter(dat, year == 2014, day == 0)
    r16 <- dplyr::filter(dat, year == 2016, day == 0)
    # One row per simulation replicate, matched by sim_i
    tibble(sim_i = r16$sim_i) |>
      mutate(
        Ratio_2014   = map_dbl(sim_i, function(s) {
          row <- r14[r14$sim_i == s, ]
          if (nrow(row) > 0 && row$ant_total[1] > 0) row$inf_total[1] / row$ant_total[1] else NA_real_
        }),
        Antigen_2016 = if_else(r16$pop_total > 0, r16$ant_total / r16$pop_total * 100, NA_real_),
        MF_2016      = if_else(r16$pop_total > 0, r16$inf_total / r16$pop_total * 100, NA_real_)
      ) |>
      dplyr::select(-sim_i)
  }, error = function(e) {
    warning("Failed to parse ", path, ": ", conditionMessage(e))
    tibble(Ratio_2014 = NA_real_, Antigen_2016 = NA_real_, MF_2016 = NA_real_)
  })
}

# ── Helper: write ABC-GLM output files matching existing format ─────────────────
write_abc_glm_files <- function(abc_fit, out_dir, label, prior) {
  write_rds(abc_fit, file.path(out_dir, "abc_fit.rds"))
  
  adj <- as.data.frame(abc_fit$adj.values)
  wts <- abc_fit$weights / sum(abc_fit$weights)

  idx  <- sample(nrow(adj), N_POSTERIOR, replace = TRUE, prob = wts)
  post <- adj[idx, ]
  colnames(post) <- c("T1", "W", "k")

  writeLines(format(post$T1, scientific = TRUE), file.path(out_dir, "Theta1.txt"))
  writeLines(format(post$k,  scientific = TRUE), file.path(out_dir, "Agg.txt"))
  writeLines(format(post$W,  scientific = TRUE), file.path(out_dir, "Work.txt"))

  theta2_val <- THETA2_VALS[THETA2_LABELS == label]
  write_csv(tibble(
    Theta_1   = median(post$T1),
    Theta_2   = theta2_val,
    Agg       = median(post$k),
    WorktoNot = median(post$W)
  ), file.path(out_dir, "TranParams.csv"))

  # Posterior characteristics
  mode_val  <- function(x) { d <- density(x); d$x[which.max(d$y)] }
  post_mat  <- as.matrix(post)
  hpd50 <- HDInterval::hdi(post_mat, credMass = 0.50)
  hpd90 <- HDInterval::hdi(post_mat, credMass = 0.90)
  hpd95 <- HDInterval::hdi(post_mat, credMass = 0.95)
  hpd99 <- HDInterval::hdi(post_mat, credMass = 0.99)

  chars <- bind_rows(
    tibble(what = "mode",                    T1 = mode_val(post$T1),        W = mode_val(post$W),        k = mode_val(post$k)),
    tibble(what = "mean",                    T1 = mean(post$T1),            W = mean(post$W),            k = mean(post$k)),
    tibble(what = "median",                  T1 = median(post$T1),          W = median(post$W),          k = median(post$k)),
    tibble(what = "quantile_50_lower_bound", T1 = quantile(post$T1, 0.25),  W = quantile(post$W, 0.25),  k = quantile(post$k, 0.25)),
    tibble(what = "quantile_50_upper_bound", T1 = quantile(post$T1, 0.75),  W = quantile(post$W, 0.75),  k = quantile(post$k, 0.75)),
    tibble(what = "quantile_90_lower_bound", T1 = quantile(post$T1, 0.05),  W = quantile(post$W, 0.05),  k = quantile(post$k, 0.05)),
    tibble(what = "quantile_90_upper_bound", T1 = quantile(post$T1, 0.95),  W = quantile(post$W, 0.95),  k = quantile(post$k, 0.95)),
    tibble(what = "quantile_95_lower_bound", T1 = quantile(post$T1, 0.025), W = quantile(post$W, 0.025), k = quantile(post$k, 0.025)),
    tibble(what = "quantile_95_upper_bound", T1 = quantile(post$T1, 0.975), W = quantile(post$W, 0.975), k = quantile(post$k, 0.975)),
    tibble(what = "quantile_99_lower_bound", T1 = quantile(post$T1, 0.005), W = quantile(post$W, 0.005), k = quantile(post$k, 0.005)),
    tibble(what = "quantile_99_upper_bound", T1 = quantile(post$T1, 0.995), W = quantile(post$W, 0.995), k = quantile(post$k, 0.995)),
    tibble(what = "HPD_50_lower_bound",      T1 = hpd50["lower","T1"], W = hpd50["lower","W"], k = hpd50["lower","k"]),
    tibble(what = "HPD_50_upper_bound",      T1 = hpd50["upper","T1"], W = hpd50["upper","W"], k = hpd50["upper","k"]),
    tibble(what = "HPD_90_lower_bound",      T1 = hpd90["lower","T1"], W = hpd90["lower","W"], k = hpd90["lower","k"]),
    tibble(what = "HPD_90_upper_bound",      T1 = hpd90["upper","T1"], W = hpd90["upper","W"], k = hpd90["upper","k"]),
    tibble(what = "HPD_95_lower_bound",      T1 = hpd95["lower","T1"], W = hpd95["lower","W"], k = hpd95["lower","k"]),
    tibble(what = "HPD_95_upper_bound",      T1 = hpd95["upper","T1"], W = hpd95["upper","W"], k = hpd95["upper","k"]),
    tibble(what = "HPD_99_lower_bound",      T1 = hpd99["lower","T1"], W = hpd99["lower","W"], k = hpd99["lower","k"]),
    tibble(what = "HPD_99_upper_bound",      T1 = hpd99["upper","T1"], W = hpd99["upper","W"], k = hpd99["upper","k"])
  )
  write_tsv(chars, file.path(out_dir, "ABC_GLM_PosteriorCharacteristics_Obs0.txt"))

  # 100-point density grid
  d_t1 <- density(post$T1, n = 100)
  d_w  <- density(post$W,  n = 100)
  d_k  <- density(post$k,  n = 100)
  write_tsv(tibble(
    number    = 1:100,
    T1        = d_t1$x, T1_density = d_t1$y,
    W         = d_w$x,  W_density  = d_w$y,
    k         = d_k$x,  k_density  = d_k$y
  ), file.path(out_dir, "ABC_GLM_PosteriorEstimates_Obs0.txt"))

  # L1 distance prior vs posterior
  l1 <- function(a, b) {
    grid  <- seq(min(a, b), max(a, b), length.out = 512)
    da    <- density(a, from = grid[1], to = grid[512], n = 512)$y
    db    <- density(b, from = grid[1], to = grid[512], n = 512)$y
    mean(abs(da - db))
  }
  write_tsv(tibble(
    T1 = l1(runif(N_POSTERIOR, T1_MIN, T1_MAX), post$T1),
    W  = l1(runif(N_POSTERIOR, W_MIN,  W_MAX),  post$W),
    k  = l1(runif(N_POSTERIOR, K_MIN,  K_MAX),  post$k)
  ), file.path(out_dir, "ABC_GLM_L1DistancePriorPosterior.txt"))

  # Prior (black) / posterior (red) density plots
  tryCatch({
    pdf(file.path(out_dir, "ABC_GLM_PosteriorPlots_Obs0.pdf"))
    plot(abc_fit, param = prior, ask = FALSE, subsample = min(1000, nrow(prior)))
    dev.off()
  }, error = function(e) {
    dev.off()
    warning("Prior/posterior plot failed (too few particles?): ", conditionMessage(e))
  })
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

  particles <- tibble(
    Sim          = integer(),
    T1           = double(),
    W            = double(),
    k            = double(),
    Ratio_2014   = double(),
    Antigen_2016 = double(),
    MF_2016      = double()
  )

  for (p in seq_len(N_PARTICLES)) {

    t1 <- runif(1, T1_MIN, T1_MAX)
    k  <- runif(1, K_MIN,  K_MAX)
    w  <- runif(1, W_MIN,  W_MAX)

    write_csv(
      tibble(Theta_1 = t1, Theta_2 = theta2, Agg = k, WorktoNot = w),
      file.path(data_dir, "TranParams.csv")
    )

    id       <- sprintf("%s_%s_p%04d", TSV_PREFIX, suffix, p)
    out_file <- file.path(output_dir, id)

    run_model(id)

    if (file.exists(out_file)) {
      stats <- parse_output(out_file)
      file.remove(out_file)
      particles <- bind_rows(particles,
        mutate(stats, Sim = p, T1 = t1, W = w, k = k)
      )
    } else {
      message(sprintf("  Warning: no output for p=%d", p))
    }

    if (p %% 5 == 0) message(sprintf("  %d / %d particles", p, N_PARTICLES))
  }

  tsv_name <- sprintf("fit_%s%s.tsv", TSV_PREFIX, suffix)
  tsvr_name <- sprintf("fit_%s%sr.tsv", TSV_PREFIX, suffix)

  write_tsv(particles, file.path(out_theta, tsv_name))
  write_tsv(dplyr::select(particles, -Ratio_2014), file.path(out_theta, tsvr_name))

  # ── ABC-GLM ──────────────────────────────────────────────────────────────────
  complete <- particles |>
    drop_na(Antigen_2016, MF_2016) |>
    group_by(Sim, T1, W, k) |>
    summarise(
      Antigen_2016 = mean(Antigen_2016),
      MF_2016      = mean(MF_2016),
      .groups = "drop"
    )
  
  message(sprintf("  %d / %d complete rows for ABC", nrow(complete), nrow(particles)))
  
  # L1 rejection: keep closest ABC_TOL fraction by L1 distance
  obs_vec <- c(OBS_ANT_2016, OBS_MF_2016)
  ss_mat  <- as.matrix(dplyr::select(complete, Antigen_2016, MF_2016))
  l1_dist <- rowSums(abs(sweep(ss_mat, 2, obs_vec)))
  accepted <- complete[l1_dist <= quantile(l1_dist, 0.03), ]
  message(sprintf("  %d particles accepted after L1 rejection (tol = %.0f%%)",
                  nrow(accepted), 0.03 * 100))

  # GLM post-sampling adjustment on accepted particles
  fit <- abc(
    target  = obs_vec,
    param   = dplyr::select(accepted, T1, W, k),
    sumstat = dplyr::select(accepted, Antigen_2016, MF_2016),
    tol     = 1.0,
    method  = "loclinear",
    transf  = "none"
  )

  prior <- tibble(
    T1 = runif(N_POSTERIOR, T1_MIN, T1_MAX),
    W  = runif(N_POSTERIOR, W_MIN,  W_MAX),
    k  = runif(N_POSTERIOR, K_MIN,  K_MAX)
  )
  write_abc_glm_files(fit, out_theta, label, prior)
  message(sprintf("  Written to %s", out_theta))
}

message("\nDone. Results in ", fit_base)
