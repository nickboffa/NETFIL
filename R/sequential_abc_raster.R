library(EasyABC)    # ABC_sequential() — sequential ABC / ABC-SMC implementation

source("R/abc_common.R")

# ── ABC settings specific to sequential ABC (ABC-SMC) ───────────────────────────
# N_PARTICLES here is the size of the final retained population — unlike
# abc_raster.R's N_PARTICLES (total draws from the prior), see cost note below.
N_PARTICLES <- 50

# EasyABC::ABC_sequential() replaces abc_raster.R's one-shot rejection step
# with Lenormand et al. (2012)'s adaptive population Monte Carlo: draw from
# the prior, then repeatedly resample + perturb + re-filter against a
# shrinking tolerance until a step's acceptance rate drops below
# SMC_P_ACC_MIN. Chosen since it needs no hand-tuned tolerance schedule and
# requires uniform priors, which is all we have here.
#
# nb_simul = N_PARTICLES / SMC_ALPHA draws per step, of which SMC_ALPHA
# survive as the retained population. Total cost isn't fixed in advance — it
# can run well past a rough estimate, since with few draws/step the
# acceptance-rate estimate is coarse and can hover near SMC_P_ACC_MIN for a
# while. Raising SMC_P_ACC_MIN stops it sooner, trading away some
# final-tolerance sharpness.
SMC_ALPHA     <- 0.5   # EasyABC default: fraction of each step's draws kept
# Back down to EasyABC's own default (0.05) from the 0.15 first used here —
# the first run's p_acc plateaued at ~0.55-0.6 for several steps before
# finally starting to drop (only reaching 0.15 at step 12, stopping at step
# 17 when it hit 0.03), so 0.15 was stopping well before the population had
# actually finished tightening. Lower means more steps/cost, but the
# tightened priors mean each step is already well-targeted.
SMC_P_ACC_MIN <- 0.05

# ── Helper: build the model() function ABC_sequential() calls per particle ─────
# EasyABC's cluster mode dispatches `model` to bare worker processes with no
# access to this session's globals, so everything model() needs has to live
# in its own closure (build_model()'s arguments and nested helpers below)
# rather than the top-level script environment (including R/abc_common.R's
# icc_to_sigma_beta()/run_model()/parse_output(), all duplicated below for the
# same reason), and every tidyverse call inside must be namespaced since
# nothing is library()'d on the worker.
#
# use_seed = TRUE hands model() a per-call integer seed as x[1], used here as
# a collision-free id for run_model()'s output file.
build_model <- function(model_dir, output_dir, theta2, suffix, N_REPS, ANT_0, fit_stats,
                         fit_year, fit_icc, fixed_sigma_group, fixed_beta_0,
                         fit_t1, fixed_t1, fit_w, fixed_w, fit_k, fixed_k,
                         fit_ant_loss, fixed_ant_loss,
                         fit_p_mda, fixed_p_mda,
                         fit_ster_to_kill, fixed_ster_to_kill,
                         fit_mda_effect, fixed_mda_effect) {
  # force() every argument: left as lazy promises, they'd still be unresolved
  # references into this session's globalenv when serialized to a worker,
  # which has its own empty globalenv — forcing bakes in the actual value.
  force(model_dir); force(output_dir); force(theta2); force(suffix)
  force(N_REPS); force(ANT_0); force(fit_stats)
  force(fit_year); force(fit_icc); force(fixed_sigma_group); force(fixed_beta_0)
  force(fit_t1); force(fixed_t1); force(fit_w); force(fixed_w); force(fit_k); force(fixed_k)
  force(fit_ant_loss); force(fixed_ant_loss)
  force(fit_p_mda); force(fixed_p_mda)
  force(fit_ster_to_kill); force(fixed_ster_to_kill)
  force(fit_mda_effect); force(fixed_mda_effect)

  # Deliberately duplicated from R/icc_to_sigma_beta.R rather than called by
  # name or re-sourced here: workers get a bare/empty globalenv with no
  # access to this session's sourced files any more than its variables (same
  # reasoning as the force() calls above), so the function body itself has
  # to travel inside the closure.
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

  # Captures stderr (e.g. model/sim.cpp's seeding-attempts-exceeded message)
  # instead of discarding it, so a particle that can't seed shows up as a
  # real R warning with the actual C++ error text, not a silent rejection.
  run_model <- function(id, t1, t2, k, w, sigma_group, beta_0,
                         ant_loss = NA_real_, p_mda = NA_real_, ster_to_kill = NA_real_,
                         mda_total_effect = NA_real_) {
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(model_dir)
    # output_epidemics() appends to this file — remove any leftover from a
    # crashed prior run so its rows don't get silently mixed in.
    unlink(file.path(output_dir, id))
    args <- c(id, N_REPS, t1, t2, k, w, sigma_group, beta_0)
    if (!is.na(ant_loss)) args <- c(args, ant_loss)               # 9th CLI arg
    if (!is.na(p_mda)) args <- c(args, p_mda)                     # 10th CLI arg — needs 9th present too
    if (!is.na(ster_to_kill)) args <- c(args, ster_to_kill)       # 11th CLI arg — needs 9th and 10th present too
    if (!is.na(mda_total_effect)) args <- c(args, mda_total_effect) # 12th CLI arg — needs 9th-11th present too
    system2("./main", args = args, stdout = FALSE, stderr = TRUE)
  }

  # output_epidemics writes quarterly; we use year == fit_year, day == 0 rows.
  parse_output <- function(path) {
    tryCatch({
      # `name` is a string column; everything else is numeric
      dat <- readr::read_csv(path, col_types = readr::cols(name = "c", .default = "d"), show_col_types = FALSE)
      rY <- dplyr::filter(dat, year == fit_year, day == 0)
      if (nrow(rY) == 0) return(tibble::tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_))

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
        ANT = dplyr::if_else(rY$pop_total > 0, rY$ant_total / rY$pop_total * 100, NA_real_),
        MF  = dplyr::if_else(rY$pop_total > 0, rY$mf_total  / rY$pop_total * 100, NA_real_),
        ICC = purrr::map_dbl(seq_len(nrow(rY)), \(i) calc_icc(rY[i, ]))
      )
    }, error = function(e) {
      warning("Failed to parse ", path, ": ", conditionMessage(e))
      tibble::tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_)
    })
  }

  function(x) {
    seed <- as.integer(round(x[1]))
    idx  <- 2   # next free slot in x[], after seed
    # T1/W/k: each independently fit, or pinned to ASM's Raster550 fit
    # (fixed_t1/w/k) — see the FIT_T1/FIT_W/FIT_K comment in R/abc_common.R.
    # Order must match the prior list below (T1, W, k).
    t1 <- if (fit_t1) { v <- x[idx]; idx <- idx + 1; v } else fixed_t1
    w  <- if (fit_w)  { v <- x[idx]; idx <- idx + 1; v } else fixed_w
    k  <- if (fit_k)  { v <- x[idx]; idx <- idx + 1; v } else fixed_k
    if (fit_icc) {
      icc10 <- x[idx]; idx <- idx + 1
      cb <- icc_to_sigma_beta(icc10, ANT_0)
    } else {
      cb <- list(sigma_group = fixed_sigma_group, beta_0 = fixed_beta_0)
    }
    # ant_loss/p_mda/ster_to_kill: fit value if this run fits it, else the
    # fixed fallback baked into this closure — always a real number so
    # model/main.cpp's CLI positions stay aligned regardless of which of the
    # three is actually being fit (see FIXED_*/run_model() comments in
    # R/abc_common.R).
    ant_loss     <- if (fit_ant_loss) { v <- x[idx]; idx <- idx + 1; v } else fixed_ant_loss
    p_mda        <- if (fit_p_mda)    { v <- x[idx]; idx <- idx + 1; v } else fixed_p_mda
    ster_to_kill <- if (fit_ster_to_kill) { v <- x[idx]; idx <- idx + 1; v } else fixed_ster_to_kill
    mda_effect   <- if (fit_mda_effect)   { v <- x[idx]; idx <- idx + 1; v } else fixed_mda_effect

    id <- sprintf("%s_s%d", suffix, seed)
    stderr_lines <- run_model(id, t1, theta2, k, w, cb$sigma_group, cb$beta_0,
                               ant_loss, p_mda, ster_to_kill, mda_effect)

    out_file <- file.path(output_dir, id)
    if (!file.exists(out_file)) {
      # Model failed to produce output — most often model/sim.cpp's seeding
      # loop hit MAX_SEED_ATTEMPTS (100) because this particle's parameters
      # can't reach country.json's seeding target at all (e.g. mismatched
      # fixed clustering params). Surface the real error rather than silently
      # rejecting — return something far from any target so the particle is
      # naturally rejected rather than NA crashing the sampler.
      detail <- if (length(stderr_lines) > 0) paste(stderr_lines, collapse = " | ") else "no stderr captured"
      warning(sprintf("Particle %s produced no output (T1=%.4g W=%.4g k=%.4g): %s", id, t1, w, k, detail))
      return(rep(1e6, length(fit_stats)))
    }
    stats <- parse_output(out_file)
    file.remove(out_file)

    # A file that exists but failed to parse (e.g. corrupted/headerless from
    # a rare write race under heavy parallel load) comes back from
    # parse_output() as all-NA rows — treat it the same as the missing-file
    # case above rather than letting NA reach ABC_sequential(), which has no
    # NA handling either.
    if (all(is.na(stats$ANT)) || all(is.na(stats$MF))) {
      warning(sprintf("Particle %s produced unparseable output (T1=%.4g W=%.4g k=%.4g)", id, t1, w, k))
      return(rep(1e6, length(fit_stats)))
    }

    icc_vals <- stats$ICC
    icc_vals[is.na(icc_vals)] <- 0  # occurs when MF prevalence is really low
    result <- c(ANT = mean(stats$ANT), MF = mean(stats$MF), ICC = mean(icc_vals))
    as.numeric(result[fit_stats])
  }
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

  # Warm-up: pre-build the .init population cache so workers don't race to
  # create it once ABC_sequential() starts.
  message("  Warm-up run to build population cache...")
  warmup_cluster <- compute_cluster_params(if (FIT_ICC) (ICC10_MIN + ICC10_MAX) / 2 else NA_real_)
  warmup_t1 <- if (FIT_T1) max(T1_MIN, 1e-4) else FIXED_T1
  warmup_k  <- if (FIT_K) K_MIN_EFF else FIXED_K
  warmup_w  <- if (FIT_W) max(W_MIN, 1e-4) else FIXED_W
  t_warmup_start <- Sys.time()
  run_model("warmup", warmup_t1, theta2, warmup_k, warmup_w,
            warmup_cluster$sigma_group, warmup_cluster$beta_0)
  t_warmup_end <- Sys.time()
  suppressWarnings(file.remove(file.path(output_dir, "warmup")))
  warmup_sec <- as.numeric(difftime(t_warmup_end, t_warmup_start, units = "secs"))

  fit_stats <- if (FIT_ICC) c("ANT", "MF", "ICC") else c("ANT", "MF")
  obs_fit   <- if (FIT_ICC) c(OBS_ANT, OBS_MF, OBS_ICC) else
               c(OBS_ANT, OBS_MF)
  obs_vec   <- c(OBS_ANT, OBS_MF, OBS_ICC)  # always 3-element for visualisation; ICC is NA when not fit

  # ── Sequential ABC ────────────────────────────────────────────────────────────
  prior <- c(
    if (FIT_T1) list(c("unif", T1_MIN, T1_MAX)) else NULL,
    if (FIT_W)  list(c("unif", W_MIN,  W_MAX)) else NULL,
    if (FIT_K)  list(c("unif", K_MIN_EFF, K_MAX)) else NULL,
    if (FIT_ICC) list(c("unif", ICC10_MIN, ICC10_MAX)) else NULL,
    if (FIT_ANT_LOSS) list(c("unif", ANT_LOSS_MIN, ANT_LOSS_MAX)) else NULL,
    if (FIT_P_MDA) list(c("unif", P_MDA_MIN, P_MDA_MAX)) else NULL,
    if (FIT_STER_TO_KILL) list(c("unif", STER_TO_KILL_MIN, STER_TO_KILL_MAX)) else NULL,
    if (FIT_MDA_EFFECT) list(c("unif", MDA_EFFECT_MIN, MDA_EFFECT_MAX)) else NULL
  )
  model <- build_model(model_dir, output_dir, theta2, suffix, N_REPS, ANT_0, fit_stats,
                        FIT_YEAR, FIT_ICC, FIXED_SIGMA_GROUP, FIXED_BETA_0,
                        FIT_T1, FIXED_T1, FIT_W, FIXED_W, FIT_K, FIXED_K,
                        FIT_ANT_LOSS, FIXED_ANT_LOSS,
                        FIT_P_MDA, FIXED_P_MDA,
                        FIT_STER_TO_KILL, FIXED_STER_TO_KILL,
                        FIT_MDA_EFFECT, FIXED_MDA_EFFECT)

  nb_simul_arg <- ceiling(N_PARTICLES / SMC_ALPHA)
  message(sprintf("  Running ABC_sequential (target %d particles/step, alpha = %.2f, p_acc_min = %.2f, %d workers)...",
                  N_PARTICLES, SMC_ALPHA, SMC_P_ACC_MIN, N_WORKERS))
  t_smc_start <- Sys.time()
  # verbose=TRUE writes EasyABC's per-step diagnostic files into the current
  # working directory at call time — setwd() into out_theta first so they
  # land there instead of littering the project root.
  old_wd <- setwd(out_theta)
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
    verbose             = TRUE
  )
  setwd(old_wd)
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
  # The final population already passed the algorithm's own tolerance, so
  # (unlike abc_raster.R's separate draw/reject split) "all particles" and
  # "accepted" below are the same set — passed twice to keep
  # write_abc_glm_files() unchanged.
  colnames(smc$param) <- PARAM_NAMES
  colnames(smc$stats)  <- fit_stats

  complete <- as_tibble(smc$param) |>
    bind_cols(as_tibble(smc$stats)) |>
    mutate(Sim = row_number(), weight = smc$weights)

  ss_mat  <- as.matrix(dplyr::select(complete, all_of(fit_stats)))
  l1_dist <- rowSums(abs(sweep(ss_mat, 2, obs_fit)))

  message(sprintf("  %d particles in final population", nrow(complete)))

  # ── GLM post-sampling adjustment ───────────────────────────────────────────────
  # Resample the weighted population down to an unweighted set (particles
  # with more SMC weight appear more often) so it drops into the same
  # local-linear regression-adjustment step abc_raster.R uses.
  resample_idx <- sample(nrow(complete), nrow(complete), replace = TRUE, prob = complete$weight)
  accepted     <- complete[resample_idx, ]

  fit <- abc(
    target  = obs_fit,
    param   = dplyr::select(accepted, all_of(PARAM_NAMES)),
    sumstat = dplyr::select(accepted, all_of(fit_stats)),
    tol     = 1.0,
    method  = "loclinear",
    transf  = "none"
  )

  prior_draws <- tibble(.rows = N_POSTERIOR)
  if (FIT_T1) prior_draws$T1 <- runif(N_POSTERIOR, T1_MIN, T1_MAX)
  if (FIT_W)  prior_draws$W  <- runif(N_POSTERIOR, W_MIN,  W_MAX)
  if (FIT_K)  prior_draws$k  <- runif(N_POSTERIOR, K_MIN,  K_MAX)
  if (FIT_ICC) prior_draws$ICC10 <- runif(N_POSTERIOR, ICC10_MIN, ICC10_MAX)
  if (FIT_ANT_LOSS) prior_draws$AntLoss <- runif(N_POSTERIOR, ANT_LOSS_MIN, ANT_LOSS_MAX)
  if (FIT_P_MDA) prior_draws$PMda <- runif(N_POSTERIOR, P_MDA_MIN, P_MDA_MAX)
  if (FIT_STER_TO_KILL) prior_draws$SterToKill <- runif(N_POSTERIOR, STER_TO_KILL_MIN, STER_TO_KILL_MAX)
  if (FIT_MDA_EFFECT) prior_draws$MdaEffect <- runif(N_POSTERIOR, MDA_EFFECT_MIN, MDA_EFFECT_MAX)

  write_csv(complete, csv_path)
  write_abc_glm_files(fit, out_theta, label, prior_draws,
                      complete, accepted, l1_dist, obs_vec)

  message(sprintf("  Written to %s", out_theta))
}

message("\nDone. Results in ", fit_base)
