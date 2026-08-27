library(furrr)
library(progressr)

handlers(global = TRUE)
handlers("txtprogressbar")  # live "N/N particles done" bar during future_map()

source("R/abc_common.R")

# ── ABC settings specific to plain rejection ABC ────────────────────────────────
# N_PARTICLES here is the total number of draws from the prior — unlike
# sequential_abc_raster.R's N_PARTICLES (the size of the final retained
# population), most of these get discarded by the ABC_TOL L1 cutoff below.
N_PARTICLES <- 60 # 600     # was 500

# ── Main: loop over theta2 values ─────────────────────────────────────────────
for (i in seq_along(THETA2_VALS)) {

  theta2 <- THETA2_VALS[i]
  label  <- THETA2_LABELS[i]
  suffix <- THETA2_SUFFIX[i]
  out_theta <- file.path(fit_base, label)
  dir.create(out_theta, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n── %s  %s (theta2 = %s) — rejection ABC ──────────────────────────",
                  SCALE, label, theta2))

  csv_name <- sprintf("fit_%s.csv", suffix)
  csv_path <- file.path(out_theta, csv_name)

  # Stack onto a previous run's particles if fit_<suffix>.csv already exists,
  # so repeated runs build up one larger pool instead of overwriting it. New
  # particles' Sim ids continue after the existing max so group_by(Sim, ...)
  # below doesn't collide old and new particles together.
  previous_particles <- NULL
  sim_offset <- 0
  if (file.exists(csv_path)) {
    previous_particles <- read_csv(csv_path, show_col_types = FALSE)
    missing_cols <- setdiff(c("Sim", "ANT", "MF", "ICC", PARAM_NAMES), names(previous_particles))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "%s is missing column(s) %s — likely from a run with different FIT_* settings. Move or delete it before stacking a new run on top.",
        csv_path, paste(missing_cols, collapse = ", ")
      ))
    }
    # Any parameter that stops being ABC-fit particle-by-particle (its own
    # FIT_* is FALSE, pinned to a FIXED_* value instead) risks pooling two
    # different experiments into one GLM regression if an old file's column
    # for it was freely fit under a different config. Old rows still carry
    # the column (just not required by PARAM_NAMES above when not fit), so
    # check each fixed one actually matches rather than silently mixing.
    fixed_param_checks <- list(
      T1         = list(fit = FIT_T1,            fixed = FIXED_T1),
      W          = list(fit = FIT_W,             fixed = FIXED_W),
      k          = list(fit = FIT_K,             fixed = FIXED_K),
      AntLoss    = list(fit = FIT_ANT_LOSS,       fixed = FIXED_ANT_LOSS),
      PMda       = list(fit = FIT_P_MDA,          fixed = FIXED_P_MDA),
      SterToKill = list(fit = FIT_STER_TO_KILL,   fixed = FIXED_STER_TO_KILL),
      MdaEffect  = list(fit = FIT_MDA_EFFECT,      fixed = FIXED_MDA_EFFECT)
    )
    for (pname in names(fixed_param_checks)) {
      spec <- fixed_param_checks[[pname]]
      if (!spec$fit && pname %in% names(previous_particles)) {
        if (!isTRUE(all.equal(unique(previous_particles[[pname]]), spec$fixed))) {
          stop(sprintf(
            "%s has %s values that don't match the current FIXED_%s (%.4g) — it's from a run where %s was freely fit. Move or delete it before stacking a FIT_%s=FALSE run on top.",
            csv_path, pname, pname, spec$fixed, pname, pname
          ))
        }
      }
    }
    # AntLoss changed units (raw daily retention probability, always < 1) ->
    # antigen half-life in days (always > 1, in practice tens to hundreds) —
    # a file from before that change would silently pool incompatible units
    # into the same regression column, since the checks above only catch a
    # FIXED_* mismatch, not a unit mismatch on a column that's still fit.
    if (FIT_ANT_LOSS && "AntLoss" %in% names(previous_particles) &&
        any(previous_particles$AntLoss < 1, na.rm = TRUE)) {
      stop(sprintf(
        "%s has AntLoss values below 1 — those are old-units (daily retention probability), not the current half-life-in-days parametrisation. Move or delete it before stacking a new run on top.",
        csv_path
      ))
    }
    sim_offset <- max(previous_particles$Sim)
    message(sprintf("  Found %d rows from previous run(s) in %s — stacking on top (Sim offset %d)",
                    nrow(previous_particles), csv_name, sim_offset))
  }

  # Warm-up: pre-build the .init population cache so parallel workers don't
  # race to create it.
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

  plan(multisession, workers = N_WORKERS)
  message(sprintf("  Running %d particles in parallel (%d workers)...", N_PARTICLES, N_WORKERS))
  t_par_start <- Sys.time()

  results_list <- future_map(seq_len(N_PARTICLES), function(p) {
    # T1/W/k: each independently fit, or pinned to ASM's Raster550 fit
    # (FIXED_T1/W/K) — see the FIT_T1/FIT_W/FIT_K comment in R/abc_common.R.
    t1       <- if (FIT_T1) runif(1, T1_MIN, T1_MAX) else FIXED_T1
    k        <- if (FIT_K) runif(1, K_MIN_EFF, K_MAX) else FIXED_K
    w        <- if (FIT_W) runif(1, W_MIN, W_MAX) else FIXED_W
    icc10    <- if (FIT_ICC) runif(1, ICC10_MIN, ICC10_MAX) else NA_real_
    # ant_loss/p_mda: fit value if this country fits them, else the fixed
    # fallback — always a real number so model/main.cpp's CLI positions stay
    # aligned regardless of which of the three is actually being fit (see
    # FIXED_*/run_model() comments in R/abc_common.R).
    ant_loss     <- if (FIT_ANT_LOSS) runif(1, ANT_LOSS_MIN, ANT_LOSS_MAX) else FIXED_ANT_LOSS
    p_mda        <- if (FIT_P_MDA) runif(1, P_MDA_MIN, P_MDA_MAX) else FIXED_P_MDA
    ster_to_kill <- if (FIT_STER_TO_KILL) runif(1, STER_TO_KILL_MIN, STER_TO_KILL_MAX) else FIXED_STER_TO_KILL
    mda_effect   <- if (FIT_MDA_EFFECT) runif(1, MDA_EFFECT_MIN, MDA_EFFECT_MAX) else FIXED_MDA_EFFECT
    cb       <- compute_cluster_params(icc10)

    id       <- sprintf("%s_p%04d", suffix, p)
    out_file <- file.path(output_dir, id)

    stderr_lines <- run_model(id, t1, theta2, k, w, cb$sigma_group, cb$beta_0,
                               ant_loss, p_mda, ster_to_kill, mda_effect)

    if (file.exists(out_file)) {
      stats <- parse_output(out_file)
      file.remove(out_file)
      row <- dplyr::mutate(stats, Sim = p + sim_offset, T1 = t1, W = w, k = k)
      if (FIT_ICC) row$ICC10 <- icc10
      if (FIT_ANT_LOSS) row$AntLoss <- ant_loss
      if (FIT_P_MDA) row$PMda <- p_mda
      if (FIT_STER_TO_KILL) row$SterToKill <- ster_to_kill
      if (FIT_MDA_EFFECT) row$MdaEffect <- mda_effect
      row
    } else {
      # Model failed to produce output — most often model/sim.cpp's seeding
      # loop hit MAX_SEED_ATTEMPTS (100) because this particle's parameters
      # can't reach country.json's seeding target at all. Surface the real
      # error rather than silently dropping the row.
      detail <- if (length(stderr_lines) > 0) paste(stderr_lines, collapse = " | ") else "no stderr captured"
      warning(sprintf("Particle %s produced no output (T1=%.4g W=%.4g k=%.4g): %s", id, t1, w, k, detail))
      NULL
    }
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)

  plan(sequential)
  t_par_end <- Sys.time()

  wall_sec    <- as.numeric(difftime(t_par_end,   t_par_start,   units = "secs"))
  warmup_sec  <- as.numeric(difftime(t_warmup_end, t_warmup_start, units = "secs"))
  compute_sec <- wall_sec * N_WORKERS   # total core-seconds across all workers

  # Append rather than overwrite, same as fit_<suffix>.csv above, so
  # timing.csv keeps one row per run instead of just the latest.
  timing_path <- file.path(out_theta, "timing.csv")
  timing_row <- tibble(
    scale        = SCALE,
    theta2       = theta2,
    n_particles  = N_PARTICLES,
    n_reps       = N_REPS,
    n_workers    = N_WORKERS,
    warmup_sec   = round(warmup_sec,  1),
    wall_sec     = round(wall_sec,    1),
    compute_sec  = round(compute_sec, 1),
    started_at   = format(t_par_start, "%Y-%m-%d %H:%M:%S"),
    finished_at  = format(t_par_end,   "%Y-%m-%d %H:%M:%S")
  )
  if (file.exists(timing_path)) {
    # started_at/finished_at are formatted timestamps, not numbers — pin them
    # to character or read_csv() guesses <datetime> and bind_rows() then
    # can't combine them with the character values in this run's new row.
    old_timing <- read_csv(
      timing_path,
      col_types = cols(started_at = "c", finished_at = "c", .default = col_guess())
    )
    timing_row <- bind_rows(old_timing, timing_row)
  }
  write_csv(timing_row, timing_path)
  message(sprintf("  Wall time: %.0fs | Compute time: %.0fs (%d workers)",
                  wall_sec, compute_sec, N_WORKERS))

  new_particles <- bind_rows(results_list)
  particles <- bind_rows(previous_particles, new_particles)
  write_csv(particles, csv_path)
  message(sprintf("  %d total rows in %s (%d new)", nrow(particles), csv_name, nrow(new_particles)))

  # ── ABC-GLM ──────────────────────────────────────────────────────────────────
  fit_stats <- if (FIT_ICC) c("ANT", "MF", "ICC") else c("ANT", "MF")
  obs_fit   <- if (FIT_ICC) c(OBS_ANT, OBS_MF, OBS_ICC) else
               c(OBS_ANT, OBS_MF)
  obs_vec   <- c(OBS_ANT, OBS_MF, OBS_ICC)  # always 3-element for visualisation; ICC is NA when not fit

  complete <- particles |>
    mutate(ICC = replace_na(ICC, 0)) |> # occurs when MF prevalence is really low
    group_by(Sim, across(all_of(PARAM_NAMES))) |>
    summarise(
      ANT = mean(ANT),
      MF  = mean(MF),
      ICC = mean(ICC, na.rm = TRUE),
      .groups = "drop"
    )

  # A rep whose output file exists but failed to parse (e.g. a corrupted/
  # headerless write under heavy parallel load — rare, but the quantile()
  # call below has no NA handling and would otherwise kill the whole run
  # after all the compute above already finished) comes back from
  # parse_output() as NA rather than missing entirely. Drop it here exactly
  # like a genuinely missing output file already is.
  n_before_na_drop <- nrow(complete)
  complete <- dplyr::filter(complete, !is.na(ANT), !is.na(MF))
  n_dropped <- n_before_na_drop - nrow(complete)
  if (n_dropped > 0) {
    warning(sprintf(
      "Dropped %d particle(s) with unparseable output (corrupted/incomplete output file). Continuing with %d.",
      n_dropped, nrow(complete)
    ))
  }

  message(sprintf("  %d / %d complete rows for ABC", nrow(complete), nrow(particles)))

  # L1 rejection: keep closest ABC_TOL fraction by L1 distance
  ss_mat  <- as.matrix(dplyr::select(complete, all_of(fit_stats)))
  l1_dist <- rowSums(abs(sweep(ss_mat, 2, obs_fit)))
  accepted <- complete[l1_dist <= quantile(l1_dist, ABC_TOL), ]
  message(sprintf("  %d particles accepted after L1 rejection (tol = %.0f%%)",
                  nrow(accepted), ABC_TOL * 100))

  # GLM post-sampling adjustment on accepted particles
  fit <- abc(
    target  = obs_fit,
    param   = dplyr::select(accepted, all_of(PARAM_NAMES)),
    sumstat = dplyr::select(accepted, all_of(fit_stats)),
    tol     = 1.0,
    method  = "loclinear",
    transf  = "none"
  )

  prior <- tibble(.rows = N_POSTERIOR)
  if (FIT_T1) prior$T1 <- runif(N_POSTERIOR, T1_MIN, T1_MAX)
  if (FIT_W)  prior$W  <- runif(N_POSTERIOR, W_MIN,  W_MAX)
  if (FIT_K)  prior$k  <- runif(N_POSTERIOR, K_MIN,  K_MAX)
  if (FIT_ICC) prior$ICC10 <- runif(N_POSTERIOR, ICC10_MIN, ICC10_MAX)
  if (FIT_ANT_LOSS) prior$AntLoss <- runif(N_POSTERIOR, ANT_LOSS_MIN, ANT_LOSS_MAX)
  if (FIT_P_MDA) prior$PMda <- runif(N_POSTERIOR, P_MDA_MIN, P_MDA_MAX)
  if (FIT_STER_TO_KILL) prior$SterToKill <- runif(N_POSTERIOR, STER_TO_KILL_MIN, STER_TO_KILL_MAX)
  if (FIT_MDA_EFFECT) prior$MdaEffect <- runif(N_POSTERIOR, MDA_EFFECT_MIN, MDA_EFFECT_MAX)

  write_abc_glm_files(fit, out_theta, label, prior,
                      complete, accepted, l1_dist, obs_vec)

  message(sprintf("  Written to %s", out_theta))
}

message("\nDone. Results in ", fit_base)
