library(tidyverse)
library(abc)        # GLM regression-adjustment step, shared by both ABC methods
library(HDInterval)
library(jsonlite)

# Shared setup for R/abc_raster.R (plain rejection ABC) and
# R/sequential_abc_raster.R (ABC-SMC via EasyABC::ABC_sequential). Each driver
# sources this file, then defines its own N_PARTICLES and runs its own
# sampling stage — that's the one part that differs enough between the two
# methods (single draw-and-reject vs. adaptive population Monte Carlo) that
# it isn't shared here.

# ── Country / scale configuration ─────────────────────────────────────────────
# COUNTRY_CONFIG below drives all country-specific settings.
COUNTRY <- "WSM"
SCALE   <- "Synth550_v1"

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
# A stale ABC.init from a previous scale would crash (wrong group_blocks
# loaded). clean_inputs.sh's paths are relative to model/, so it must run
# with model/ as the working directory.
message("Clearing stale population cache...")
old_wd <- setwd(model_dir)
system2("bash", args = "clean_inputs.sh", stdout = FALSE, stderr = FALSE)
setwd(old_wd)

# ── Per-country fit configuration ───────────────────────────────────────────────
# FIT_YEAR must be <= country.json's end_year, or parse_output() below never
# finds a matching row and every particle comes back NA. WSM's 2018 survey is
# reserved for validating the fit afterwards, not for fitting to here.
#
# FIT_ICC: whether ICC10 (the 2010 seeding ICC) is a free ABC parameter, fit
# against OBS_ICC. WSM has no ICC data, so its clustering params
# (sigma_group, beta_0) are instead held fixed from
# data/WSM/Scales/<SCALE>/clustering_params.csv.
COUNTRY_CONFIG <- list(
  ASM = list(
    FIT_YEAR  = 2016,
    OBS_ANT   = 6.2,   # ANT prevalence % (Lau et al. 2020, community survey age ≥8)
    OBS_MF    = 1.59,  # MF prevalence % (25.6% of ANT+ Mf-positive, back-calculated)
    OBS_ICC   = 0.49,  # 0.5134347 # exponential y = a * exp (b * x)
    FIT_ICC   = TRUE,
    ICC10_MIN = 0.4, ICC10_MAX = 0.9,
    # This is the ASM Raster550 fit WSM's FIXED_T1/K/W below borrow from —
    # ASM always fits all three itself.
    FIT_T1 = TRUE, FIT_W = TRUE, FIT_K = TRUE,
    FIT_ANT_LOSS = FALSE,
    FIT_P_MDA = FALSE,
    FIT_STER_TO_KILL = FALSE,
    FIT_MDA_EFFECT = FALSE
  ),
  WSM = list(
    FIT_YEAR  = 2004,
    OBS_ANT   = 1.13,  # Ag prevalence % (2004 survey)
    OBS_MF    = 0.39,  # Mf prevalence % (2004 survey)
    OBS_ICC   = NA_real_,
    FIT_ICC   = FALSE,
    ICC10_MIN = NA_real_, ICC10_MAX = NA_real_,
    # T1 was previously fixed here: letting it chase WSM's low target ANT in
    # isolation (with k pinned high) is what drove it to the prior floor,
    # where mate-limitation collapsed MF regardless of any other parameter.
    # A manual sweep (T1 x k x MdaEffect, holding AntLoss/SterToKill/PMda at
    # their found-favourable values) found ANT=1.23%/MF=0.32%/ratio=3.8
    # (target 1.13%/0.39%/2.9) at T1=0.0015, k=0.01, MdaEffect=0.48 — T1 and k
    # moving DOWN TOGETHER is what avoids the collapse, since low k lets a
    # much lower population-level transmission intensity still produce
    # enough mated pairs. T1 is re-enabled now that T1_MIN/MAX (below) are
    # tightened around that jointly-found region rather than the old
    # T1-in-isolation range. W stays fixed — the sweep never found it to
    # matter and never needed to move it from ASM's fit.
    FIT_T1 = TRUE, FIT_W = FALSE, FIT_K = TRUE,
    # AntLoss was previously fit here but is parked for now while
    # ster_to_kill (below) is tried as the explanation for the low-Mf/high-Ag
    # mismatch instead — see model/mda.h.
    FIT_ANT_LOSS = TRUE,
    FIT_P_MDA = FALSE,
    # ster_to_kill: re-splits DA/IDA's combined kill+sterilise effect between
    # the two (model/mda.h's DA_IDA_TOTAL_EFFECT). Fixed at 0 (pure kill, no
    # sterilisation) rather than fit — every fit so far found low values
    # favoured (best particles clustered at 0.003-0.28, corr with fit error
    # only 0.07 i.e. barely mattered), and fixing it drops a dimension from
    # an already particle-starved GLM regression (only ~5% of particles
    # survive L1 rejection, shared across all the other free parameters).
    FIT_STER_TO_KILL = FALSE,
    # mda_total_effect: the combined probability that MDA kills OR sterilises
    # a DA/IDA worm at all (literature placeholder 0.649) — distinct from
    # SterToKill, which only re-splits this total. Lowering it means MORE
    # worms escape MDA fully fertile, which should raise Mf directly (more
    # mated pairs survive treatment intact) rather than just reshuffling
    # which "affected" bucket a worm falls into.
    FIT_MDA_EFFECT = TRUE
  )
)[[COUNTRY]]

FIT_YEAR  <- COUNTRY_CONFIG$FIT_YEAR
OBS_ANT   <- COUNTRY_CONFIG$OBS_ANT
OBS_MF    <- COUNTRY_CONFIG$OBS_MF
OBS_ICC   <- COUNTRY_CONFIG$OBS_ICC
FIT_ICC   <- COUNTRY_CONFIG$FIT_ICC
ICC10_MIN <- COUNTRY_CONFIG$ICC10_MIN
ICC10_MAX <- COUNTRY_CONFIG$ICC10_MAX
FIT_T1 <- COUNTRY_CONFIG$FIT_T1
FIT_W  <- COUNTRY_CONFIG$FIT_W
FIT_K  <- COUNTRY_CONFIG$FIT_K
FIT_ANT_LOSS <- COUNTRY_CONFIG$FIT_ANT_LOSS
FIT_P_MDA <- COUNTRY_CONFIG$FIT_P_MDA
FIT_STER_TO_KILL <- COUNTRY_CONFIG$FIT_STER_TO_KILL
FIT_MDA_EFFECT <- COUNTRY_CONFIG$FIT_MDA_EFFECT

# Prior bounds + fixed fallback for AntLoss/p_mda/ster_to_kill. Deliberately
# NOT part of COUNTRY_CONFIG above: these are prior guesses about a
# drug/coverage effect, not data specific to one country, so toggling
# FIT_ANT_LOSS/FIT_P_MDA/FIT_STER_TO_KILL is always safe on its own — the
# bounds already exist regardless of which flags are on. (Contrast with
# ICC10_MIN/MAX above, which genuinely can't be filled in for WSM without
# also supplying an OBS_ICC target — WSM has no ICC survey data at all.)
#
# Second round of tightening, this time from an actual ABC run's 240
# particles under the first tightened priors (not just the manual sweep) —
# ranked by L1 distance to target and checked for corr(param, L1 distance)
# across all 240 rows:
#   T1:      corr +0.68 (strong) — lower is better; top-40 clustered 0.0008-0.0026
#   MdaEffect: corr -0.62 (strong) — HIGHER is better here, the opposite of what
#              the earlier manual-sweep-based floor (0.3) assumed; top-40
#              clustered 0.38-0.60
#   k:       corr -0.18 (weak) — top-40 spanned almost the whole prior, left alone
#   AntLoss: corr -0.01 (none) — top-40 spanned almost the whole prior, left alone
#   SterToKill: corr +0.07 (weak, low favoured) — fixed at 0 instead (see
#               FIT_STER_TO_KILL above) rather than tightened further
# Widened per explicit user request after reading posterior.png as wanting
# ~0-70. Note: the raw simulated SMC particles (fit_1.csv, weighted p10=14.8/
# p90=36.9) never actually pressed against the old [10,45] bounds the way
# T1/k/MdaEffect did — so this isn't evidence-driven tightening/loosening the
# way those were. But a bounded prior can never show ceiling-pressure for a
# true optimum that sits outside it in the first place, so widening here is a
# reasonable hedge even without particles proving it. Floor kept just above 0
# (not exactly 0) since AntLoss feeds daily_retention_from_halflife(halflife)
# = 0.5^(1/halflife), which divides by halflife.
ANT_LOSS_MIN <- 1;  ANT_LOSS_MAX <- 70
# PMda: never varied off its ceiling (1.0, full documented historical
# coverage) in the manual sweep that first set this range — every attempt at
# lowering it only pushes prevalence further from a low target.
P_MDA_MIN    <- 0.8;  P_MDA_MAX    <- 1.0
# SterToKill's bounds stay defined (same reasoning as ANT_LOSS_MIN/MAX etc.
# above) even though FIT_STER_TO_KILL is FALSE, so re-enabling it later
# doesn't need these rediscovered.
STER_TO_KILL_MIN <- 0.0; STER_TO_KILL_MAX <- 0.3
# MdaEffect: raised the floor from 0.3 to 0.4 — the 240-particle data showed
# that floor was actively counted against fit quality (best particles never
# went below ~0.38), the opposite of the manual sweep's single-point guess.
# Ceiling raised 0.6 -> 0.7: the first sequential-ABC run's weighted p90 sat
# at 0.595, right against the old ceiling — genuine boundary pressure (all
# 30 raw simulated particles stayed within [0.4,0.6], so this isn't a GLM-
# regression-adjustment artifact, unlike the wider excursions posterior.png
# showed for AntLoss).
MDA_EFFECT_MIN   <- 0.55; MDA_EFFECT_MAX   <- 0.7

# model/main.cpp's CLI is positional: AntLoss (arg 9), p_mda (arg 10),
# ster_to_kill (arg 11), mda_total_effect (arg 12) each require every slot
# before them to be present too. Rather than track which FIT_* combination
# keeps that contiguous, the particle loops below always send all four —
# using the ABC-sampled value when a parameter is being fit, or one of the
# FIXED_* constants below otherwise — the same pattern already used for
# sigma_group/beta_0.
default_params_json <- fromJSON(file.path(data_dir, "default_params.json"))
FIXED_ANT_LOSS <- log(0.5) / log(default_params_json$daily_prob_lose_ant)  # load_config()'s own default (≈90-day half-life), echoed back explicitly
FIXED_P_MDA <- 1.0                                      # network.h default: no coverage discount
FIXED_STER_TO_KILL <- 0                                        # pure kill, no sterilisation — see FIT_STER_TO_KILL comment above
FIXED_MDA_EFFECT <- 0.649                                     # model/mda.h's DA_IDA_TOTAL_EFFECT literature placeholder

# FIXED_T1/K/W: ASM's own Raster550 fit, read fresh each run (not hardcoded)
# so it stays in sync if ASM gets refit. Only meaningful when FIT_T1/FIT_W/
# FIT_K is FALSE somewhere (currently WSM's FIT_T1/FIT_W) — read
# unconditionally since it's cheap and the file already exists regardless of
# which country is active.
asm_raster550_fit <- read_csv(
  file.path(data_dir, "ASM", "Fitted", "Raster550", "Theta_1", "tran_params.csv"),
  show_col_types = FALSE
)
FIXED_T1 <- asm_raster550_fit$Theta_1[1]
FIXED_K <- asm_raster550_fit$Agg[1]
FIXED_W <- asm_raster550_fit$WorktoNot[1]

stopifnot(
  is.finite(ANT_LOSS_MIN), is.finite(ANT_LOSS_MAX),
  is.finite(P_MDA_MIN), is.finite(P_MDA_MAX),
  is.finite(STER_TO_KILL_MIN), is.finite(STER_TO_KILL_MAX),
  is.finite(MDA_EFFECT_MIN), is.finite(MDA_EFFECT_MAX),
  is.finite(FIXED_T1), is.finite(FIXED_K), is.finite(FIXED_W)
)

# ── Semi-informative uniform priors ─────────────────────────────────────────────
# T1 tightened a second time, from the same 240-particle ABC run referenced
# in the AntLoss/PMda/SterToKill/MdaEffect comment above: corr(T1, L1
# distance) = +0.68 across all 240 rows (strong — lower is consistently
# better here), with the top-40 particles clustered 0.0008-0.0026 and none
# above ~0.0026. The old floor (0.0008) sat right at that cluster's edge, so
# only the ceiling needed to come down (0.004 -> 0.0025); a little headroom
# is kept below the observed floor in case the true optimum sits lower still.
# k left alone — corr only -0.18 (weak) and the top-40 spanned nearly this
# entire range, so there's no real evidence to narrow it further yet.
# T1_MAX raised 0.0025 -> 0.003: the first sequential-ABC run's weighted p90
# sat at 0.00249, right against the old ceiling — genuine boundary pressure
# (all 30 raw simulated particles stayed within [0.0006,0.0025], so this
# isn't a GLM-regression-adjustment artifact).
T1_MIN <- 0.002; T1_MAX <- 0.0035
K_MIN  <- 0.005;   K_MAX  <- 0.02 # don't make min too small, otherwise model will crash sometimes
W_MIN  <- 0.0;    W_MAX  <- 0.8   # untouched — W was never varied in the sweep, ASM's fixed value was never implicated

# k -> 0 makes bite_gamma(0, Inf) hang forever during agent construction, so
# the prior actually sampled from below floors K_MIN at this value.

# Free ABC parameters for this country — downstream code selects columns via
# PARAM_NAMES, so T1/W/k/ICC10/AntLoss/PMda/SterToKill drop out cleanly when
# their FIT_* flag is FALSE.
PARAM_NAMES <- c(if (FIT_T1) "T1",
                 if (FIT_W) "W",
                 if (FIT_K) "k",
                 if (FIT_ICC) "ICC10",
                 if (FIT_ANT_LOSS) "AntLoss",
                 if (FIT_P_MDA) "PMda",
                 if (FIT_STER_TO_KILL) "SterToKill",
                 if (FIT_MDA_EFFECT) "MdaEffect")

# ── ABC settings shared by both methods — each driver sets its own N_PARTICLES,
# since the two methods use it differently (total draws for rejection vs.
# final retained population size for sequential) ────────────────────────────────
N_REPS      <- 2     # was 5
N_POSTERIOR <- 10000   # was 10000, only affects output resolution, not fitted parameters
ABC_TOL     <- 0.05  # rejection: actual L1 acceptance quantile. sequential: cosmetic plot reference line only.
N_WORKERS   <- 6     # parallel workers (= cores to use; leave ≥1 free for OS/VSCode)

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

# Target 2010 seeded ANT prevalence, used by icc_to_sigma_beta() below.
country_cfg <- fromJSON(file.path(data_dir, COUNTRY, "country.json"))
ANT_0 <- country_cfg$seeding$init_ant_prev

# ── Fixed clustering params for countries without an ICC target ────────────────
FIXED_SIGMA_GROUP <- NA_real_
FIXED_BETA_0      <- NA_real_
if (!FIT_ICC) {
  fixed_cluster <- read_csv(
    file.path(data_dir, COUNTRY, "Scales", SCALE, "clustering_params.csv"),
    show_col_types = FALSE
  )
  FIXED_SIGMA_GROUP <- fixed_cluster$sigma_group[1]
  FIXED_BETA_0      <- fixed_cluster$beta_0[1]
}

# ── Helper: convert a 2010 seeding ICC into (sigma_group, beta_0) ──────────────
# See R/icc_to_sigma_beta.R. Used directly here since ICC10 is a free ABC
# parameter (no distance-decay curve involved). sequential_abc_raster.R's
# build_model() keeps its own duplicate of this inside its closure — EasyABC's
# cluster mode dispatches to bare worker processes with no access to this
# session's globals, so that copy can't just call this one.
source("R/icc_to_sigma_beta.R")

# ── Helper: cluster params for one particle — fit via ICC10, or fixed ──────────
compute_cluster_params <- function(icc10) {
  if (FIT_ICC) icc_to_sigma_beta(icc10, ANT_0)
  else list(sigma_group = FIXED_SIGMA_GROUP, beta_0 = FIXED_BETA_0)
}

# ── Helper: run model from model/ so DATADIR and OUTDIR resolve correctly ───────
# Captures stderr (e.g. model/sim.cpp's seeding-attempts-exceeded message)
# instead of discarding it, so a particle that can't seed shows up as a real
# R warning with the actual C++ error text, not a silently dropped row.
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

# ── Helper: parse output CSV for FIT_YEAR summary statistics ───────────────────
# output_epidemics writes quarterly; we use year == FIT_YEAR, day == 0 rows.
parse_output <- function(path) {
  tryCatch({
    # `name` is a string column; everything else is numeric
    dat <- read_csv(path, col_types = cols(name = "c", .default = "d"), show_col_types = FALSE)
    rY <- dplyr::filter(dat, year == FIT_YEAR, day == 0)
    if (nrow(rY) == 0) return(tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_))

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
      ANT = if_else(rY$pop_total > 0, rY$ant_total / rY$pop_total * 100, NA_real_),
      MF  = if_else(rY$pop_total > 0, rY$mf_total  / rY$pop_total * 100, NA_real_),
      ICC = purrr::map_dbl(seq_len(nrow(rY)), \(i) calc_icc(rY[i, ]))
    )
  }, error = function(e) {
    warning("Failed to parse ", path, ": ", conditionMessage(e))
    tibble(ANT = NA_real_, MF = NA_real_, ICC = NA_real_)
  })
}

# ── Helper: write ABC-GLM output files (abc_fit.rds, tran_params.csv,
# clustering_params.csv, diagnostic PNGs) ───────────────────────────────────────
write_abc_glm_files <- function(abc_fit, out_dir, label, prior,
                                 particles, accepted, l1_dist, obs_vec) {
  # ── Reload essentials ────────────────────────────────────────────────────────
  write_rds(abc_fit, file.path(out_dir, "abc_fit.rds"))

  adj <- as.data.frame(abc_fit$adj.values)
  wts <- abc_fit$weights / sum(abc_fit$weights)
  idx  <- sample(nrow(adj), N_POSTERIOR, replace = TRUE, prob = wts)
  post <- adj[idx, ]
  colnames(post) <- PARAM_NAMES

  # loclinear regression doesn't respect parameter bounds — with few accepted
  # particles it can extrapolate a posterior median outside the parameter's
  # valid range entirely (seen in practice: negative k, negative SterToKill,
  # AntLoss above its prior range). Clamp defensively rather than write out a
  # physically meaningless fitted value.
  clamp_posterior_median <- function(x, lo, hi, name) {
    if (x < lo || x > hi) {
      warning(sprintf(
        "Posterior median %s (%.4g) is outside [%.4g, %.4g] — GLM adjustment likely extrapolated past valid bounds (check accepted-particle count). Clamping to fit.",
        name, x, lo, hi
      ))
      x <- min(max(x, lo), hi)
    }
    x
  }

  theta2_val <- THETA2_VALS[THETA2_LABELS == label]
  tran_params <- tibble(
    Theta_1   = if (FIT_T1) clamp_posterior_median(median(post$T1), T1_MIN, T1_MAX, "T1") else FIXED_T1,
    Theta_2   = theta2_val,
    Agg       = if (FIT_K) clamp_posterior_median(median(post$k), K_MIN, K_MAX, "k") else FIXED_K,
    WorktoNot = if (FIT_W) clamp_posterior_median(median(post$W), W_MIN, W_MAX, "W") else FIXED_W
  )
  # kept for continuity even though clustering_params.csv is the real source
  if (FIT_ICC) tran_params$ICC10 <- clamp_posterior_median(median(post$ICC10), 0.001, 0.999, "ICC10")
  # read_parameters() (model/initialise_pop.cpp) reads columns 5-8
  # positionally as AntLoss/PMda/SterToKill/MdaEffect — once any of the four
  # is fit, all four must be written (fit median, or a FIXED_* fallback for
  # the ones this country isn't fitting) so a later column can't be misread
  # as an earlier one.
  if (FIT_ANT_LOSS || FIT_P_MDA || FIT_STER_TO_KILL || FIT_MDA_EFFECT) {
    tran_params$AntLoss    <- if (FIT_ANT_LOSS) clamp_posterior_median(median(post$AntLoss), ANT_LOSS_MIN, ANT_LOSS_MAX, "AntLoss") else FIXED_ANT_LOSS
    tran_params$PMda       <- if (FIT_P_MDA) clamp_posterior_median(median(post$PMda), P_MDA_MIN, P_MDA_MAX, "PMda") else FIXED_P_MDA
    tran_params$SterToKill <- if (FIT_STER_TO_KILL) clamp_posterior_median(median(post$SterToKill), STER_TO_KILL_MIN, STER_TO_KILL_MAX, "SterToKill") else FIXED_STER_TO_KILL
    tran_params$MdaEffect  <- if (FIT_MDA_EFFECT) clamp_posterior_median(median(post$MdaEffect), MDA_EFFECT_MIN, MDA_EFFECT_MAX, "MdaEffect") else FIXED_MDA_EFFECT
  }
  write_csv(tran_params, file.path(out_dir, "tran_params.csv"))

  if (FIT_ICC) {
    icc10_median <- tran_params$ICC10
    fitted_cluster <- icc_to_sigma_beta(icc10_median, ANT_0)
  } else {
    # No ICC target — clustering params stay fixed at what every particle used.
    fitted_cluster <- list(sigma_group = FIXED_SIGMA_GROUP, beta_0 = FIXED_BETA_0)
  }
  # icc10/ant_0 are cosmetic only — model/initialise_pop.cpp's
  # read_parameters() only ever reads the first two columns
  # (sigma_group, beta_0) off this file's data row and ignores the rest, so
  # these are purely for visually checking what sigma_group/beta_0 were
  # derived from, not consumed anywhere.
  cluster_tibble <- tibble(
    sigma_group = fitted_cluster$sigma_group,
    beta_0      = fitted_cluster$beta_0
  )
  if (FIT_ICC) {
    cluster_tibble$icc10 <- icc10_median
    cluster_tibble$ant_0 <- ANT_0
  }
  write_csv(cluster_tibble, file.path(out_dir, "clustering_params.csv"))

  # ── Diagnostic PNGs ──────────────────────────────────────────────────────────
  save_png <- function(p, name, w = 10, h = 6) {
    ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 150)
  }

  # Inter-run vs overall variation — fit_stats rather than a hardcoded
  # c(ANT, MF, ICC): sequential_abc_raster.R's `particles`/`complete` only
  # ever carries the columns in fit_stats (from smc$stats), unlike
  # abc_raster.R's, which always has ICC regardless of FIT_ICC since
  # parse_output() computes it unconditionally. Referencing fit_stats here
  # (resolved from the calling driver's global env, same as facet_stats
  # below) works for both.
  variation_dat <- particles |>
    pivot_longer(all_of(fit_stats), names_to = "stat", values_to = "value")

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

  # All particles vs accepted vs observed — one facet per fitted stat other
  # than ANT (the x-axis). ICC is only fitted for countries with FIT_ICC =
  # TRUE, so it only appears here when it was actually part of the ABC target
  # — otherwise the facet would just be empty/NA (nothing to compare).
  obs_labels  <- c(MF = "MF prevalence (%)", ICC = "ICC")
  facet_stats <- setdiff(fit_stats, "ANT")
  obs_idx     <- c(ANT = 1, MF = 2, ICC = 3)

  particle_long <- bind_rows(lapply(facet_stats, function(s)
    particles |> mutate(y = .data[[s]], facet = s)))
  accepted_long <- bind_rows(lapply(facet_stats, function(s)
    accepted |> mutate(y = .data[[s]], facet = s)))
  obs_long <- tibble(
    x     = obs_vec[1],
    y     = obs_vec[obs_idx[facet_stats]],
    facet = facet_stats
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
           x = "Antigen prevalence (%)") +
      theme_minimal() +
      theme(axis.title.y = element_blank()),
    "particles.png", w = 14
  )

  # Prior vs posterior (regression-adjusted)
  adj_df <- as.data.frame(abc_fit$adj.values) |> setNames(PARAM_NAMES)
  save_png(
    bind_rows(
      prior  |> mutate(type = "Prior"),
      adj_df |> mutate(type = "Posterior")
    ) |>
      pivot_longer(all_of(PARAM_NAMES), names_to = "param", values_to = "value") |>
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
