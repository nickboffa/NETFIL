library(tidyverse)
library(furrr)

# ── What this does ──────────────────────────────────────────────────────────
# Runs a fixed, no-MDA scenario N times (independent population + disease
# draws each time) and statistically compares the ANT/MF prevalence
# trajectory against a stored baseline. Use it after making changes to
# model/*.cpp|h to check whether behavior actually shifted, or whether any
# difference you're seeing is just run-to-run stochastic noise (the model's
# RNG is seeded from wall-clock time, so no two runs are ever identical —
# see model/rng.cpp).
#
# Usage (from repo root):
#   Rscript R/regression_test.R                    # compare current code vs baseline (n=50)
#   Rscript R/regression_test.R --n=100            # override rep count for this run
#   Rscript R/regression_test.R --fast             # reps share ONE population draw
#                                                   # (parallel, same pattern as
#                                                   # R/abc_raster.R — much faster, but
#                                                   # see WARNING below: not safe as a
#                                                   # regression check on its own)
#   Rscript R/regression_test.R --update-baseline  # regenerate tests/regression/baseline.csv
#                                                   # from the CURRENT code, at n=100 (it's
#                                                   # a one-off cost, so may as well be well-
#                                                   # powered). Do this deliberately (e.g.
#                                                   # after a real recalibration) and commit
#                                                   # the updated baseline with an explanation.
#
# Default rep counts: 100 for --update-baseline, 50 for a normal comparison
# run. n=30 was tried first and dropped — it landed WARN (d up to 0.41) on
# genuinely unchanged code purely from sampling noise, which isn't a useful
# signal. --n=N overrides whichever default applies to the mode you're running.
#
# WARNING about --fast: population construction (ages, group assignment) is
# cached and shared across all N reps within ONE script invocation — only
# disease seeding/transmission varies rep to rep. Since the baseline and the
# comparison run are always SEPARATE invocations, each draws its own single
# population realization, and that between-invocation population variance is
# NOT averaged out by adding more reps (they all share one draw). Testing this
# in practice: two --fast runs of the SAME unmodified code, n=4 each, gave
# Cohen's d up to 2.4 on ANT prevalence — a false FAIL. Default mode clears
# the population cache before every single rep (matching the methodology
# validated in this project's actual old-vs-new bisection, which needed exactly
# this to get a clean n=50 result: d=0.055-0.10, p=0.6-0.8 on truly unchanged
# code). It's slower, but it's the only mode proven not to cry wolf.
#
# Exit status: 0 if PASS/WARN, 1 if FAIL (so it can gate a merge/PR if you want).
#
# NOTE: this temporarily overwrites the "active" files in data/ (groups.csv,
# road_dist.*, country.json, tran_params.csv, clustering_params.csv,
# current_state.txt) with a fixed reference scenario, and clears the
# population cache ($config/). Your previous active config is restored and
# the cache is cleared again when the script exits (even on error).

# ── Fixed reference scenario ─────────────────────────────────────────────────
# Deliberately NOT read from the ABC-fitted files or whatever's currently
# "active" in data/ — those drift over time as fitting continues, which would
# make the baseline non-comparable across re-runs. Pinned here instead.
COUNTRY <- "ASM"
SCALE   <- "Raster550"
THETA2  <- 1
THETA2_SUFFIX <- "1"

REF_T1          <- 0.004585814591960409
REF_K           <- 0.13049485785583045
REF_W           <- 0.11835553231104337
REF_SIGMA_GROUP <- 1.1311
REF_BETA_0      <- -3.9515

# ── CLI args ─────────────────────────────────────────────────────────────────
cli_args         <- commandArgs(trailingOnly = TRUE)
UPDATE_BASELINE  <- "--update-baseline" %in% cli_args
FAST_MODE        <- "--fast" %in% cli_args   # see WARNING above — off by default on purpose
n_arg            <- grep("^--n=", cli_args, value = TRUE)
N_DEFAULT        <- if (UPDATE_BASELINE) 100 else 50  # baseline is a one-off cost — go big
N_REPS           <- if (length(n_arg) == 1) as.integer(sub("^--n=", "", n_arg)) else N_DEFAULT
N_WORKERS        <- 6

# Verdict thresholds on Cohen's d (N-independent, unlike the p-value — see the
# old-vs-new investigation this script grew out of: p=0.07 at n=20 turned into
# p=0.78 at n=50 for the SAME fixed code, purely from sampling noise).
D_WARN <- 0.3
D_FAIL <- 0.5

# ── Paths ──────────────────────────────────────────────────────────────────
script_dir  <- here::here()
model_dir   <- file.path(script_dir, "model")
data_dir    <- file.path(script_dir, "data")
output_dir  <- file.path(script_dir, "output")
config_dir  <- file.path(script_dir, "$config")
test_dir    <- file.path(script_dir, "tests", "regression")
baseline_csv  <- file.path(test_dir, "baseline.csv")
baseline_meta <- file.path(test_dir, "baseline_meta.txt")

dir.create(output_dir, showWarnings = FALSE)
dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)

# ── Back up + restore the active data/ files (this script's changes are
# temporary; the user's previous scale/config should come back untouched) ───
active_files <- c("groups.csv", "road_dist.csv", "road_dist.bin",
                   "euc_dist.csv", "euc_dist.bin", "tran_params.csv",
                   "clustering_params.csv", "country.json", "current_state.txt")
backup_dir <- file.path(tempdir(), paste0("netfil_regtest_backup_", Sys.getpid()))
dir.create(backup_dir, showWarnings = FALSE)

present_before <- active_files[file.exists(file.path(data_dir, active_files))]
if (length(present_before) > 0) {
  file.copy(file.path(data_dir, present_before), file.path(backup_dir, present_before), overwrite = TRUE)
}

restore_active_files <- function() {
  message("Restoring your previous active data/ config...")
  for (f in active_files) {
    bak <- file.path(backup_dir, f)
    live <- file.path(data_dir, f)
    if (file.exists(bak)) {
      file.copy(bak, live, overwrite = TRUE)
    } else if (file.exists(live)) {
      file.remove(live)  # script created this fresh; wasn't there before
    }
  }
  system2("bash", args = file.path(model_dir, "clean_inputs.sh"), stdout = FALSE, stderr = FALSE)
  unlink(backup_dir, recursive = TRUE)
}
on.exit(restore_active_files(), add = TRUE)

# ── Build the model fresh (this is the whole point — pick up current source) ─
message("Building model/main from current source...")
old_wd <- setwd(model_dir)
cpp_files <- Sys.glob("*.cpp")
build_ok <- system2("g++", args = c("-std=c++20", "-O2", cpp_files, "-o", "main")) == 0
setwd(old_wd)
if (!build_ok) stop("Build failed — fix compile errors before running the regression test.")

# ── Point data/ at the fixed reference scenario ──────────────────────────────
# setup.sh has no `set -e`, so a failed `cp` (e.g. a OneDrive sync hiccup —
# this repo has hit that repeatedly) prints an error but does NOT make the
# script exit non-zero, and leaves a stale file from whatever scale was
# active before silently in place. Don't trust the exit code alone — verify
# the files that actually matter (group count / topology) landed correctly,
# and retry via R's own file.copy if not.
message(sprintf("Setting up %s/%s data...", COUNTRY, SCALE))
old_wd <- setwd(model_dir)
system2("bash", args = c("setup.sh", "-country", COUNTRY, "-scale", SCALE, "-t2", THETA2_SUFFIX))
setwd(old_wd)

scale_dir <- file.path(data_dir, COUNTRY, "Scales", SCALE)
verify_and_retry_copy <- function(filename, max_tries = 5) {
  src <- file.path(scale_dir, filename)
  dst <- file.path(data_dir, filename)
  if (!file.exists(src)) return(invisible(NULL))  # e.g. road_dist.bin may not exist — fine
  for (attempt in seq_len(max_tries)) {
    if (file.exists(dst) && file.size(dst) == file.size(src)) return(invisible(NULL))
    message(sprintf("  %s missing/mismatched (attempt %d/%d) — retrying copy...",
                    filename, attempt, max_tries))
    suppressWarnings(file.remove(dst))
    file.copy(src, dst, overwrite = TRUE)
  }
  if (!file.exists(dst) || file.size(dst) != file.size(src)) {
    stop(sprintf(
      "Could not get a correct copy of %s after %d attempts (size mismatch: src=%s dst=%s). ",
      filename, max_tries,
      file.size(src), if (file.exists(dst)) file.size(dst) else "MISSING"),
      "This usually means a stale file from a DIFFERENT scale is silently in place — ",
      "do not trust results from this run.")
  }
}
for (f in c("groups.csv", "road_dist.csv", "road_dist.bin", "euc_dist.csv", "euc_dist.bin")) {
  verify_and_retry_copy(f)
}
verify_and_retry_copy_named <- function(src, dst_name, max_tries = 5) {
  dst <- file.path(data_dir, dst_name)
  for (attempt in seq_len(max_tries)) {
    if (file.exists(dst) && file.size(dst) == file.size(src)) return(invisible(NULL))
    message(sprintf("  %s missing/mismatched (attempt %d/%d) — retrying copy...",
                    dst_name, attempt, max_tries))
    suppressWarnings(file.remove(dst))
    file.copy(src, dst, overwrite = TRUE)
  }
  if (!file.exists(dst) || file.size(dst) != file.size(src)) {
    stop(sprintf("Could not get a correct copy of %s after %d attempts.", dst_name, max_tries))
  }
}
verify_and_retry_copy_named(file.path(data_dir, COUNTRY, "country.json"), "country.json")
verify_and_retry_copy_named(file.path(data_dir, COUNTRY, "initaggs.csv"), "initaggs.csv")
message("  Data files verified against source.")

# setup.sh doesn't copy clustering_params.csv — pin it to the fixed reference too.
writeLines(
  c("sigma_group,beta_0", paste(REF_SIGMA_GROUP, REF_BETA_0, sep = ",")),
  file.path(data_dir, "clustering_params.csv")
)

# ── Helper: run model from model/ so DATADIR/OUTDIR resolve correctly ───────
run_model <- function(id) {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(model_dir)
  system2("./main", args = c(id, 1, REF_T1, THETA2, REF_K, REF_W), stdout = FALSE, stderr = FALSE)
}

# ── Helper: parse one run's output CSV into a tidy ANT/MF trajectory ────────
parse_run <- function(path) {
  dat <- read_csv(path, col_types = cols(name = "c", .default = "d"), show_col_types = FALSE)
  dat |>
    filter(day %% 91 == 0) |>
    transmute(
      year, day,
      ant_prev = if_else(pop_total > 0, ant_total / pop_total * 100, NA_real_),
      mf_prev  = if_else(pop_total > 0, mf_total  / pop_total * 100, NA_real_)
    )
}

# ── Run N_REPS ────────────────────────────────────────────────────────────
if (!FAST_MODE) {
  # Default, validated mode: clear the population cache before every rep, so
  # each one is a fully independent draw (population construction included).
  # Sequential — can't safely parallelize while clearing a shared cache file.
  message(sprintf("Running %d reps sequentially (independent population each time)...", N_REPS))
  results_list <- vector("list", N_REPS)
  for (i in seq_len(N_REPS)) {
    system2("bash", args = file.path(model_dir, "clean_inputs.sh"), stdout = FALSE, stderr = FALSE)
    id       <- sprintf("regtest_%04d", i)
    out_file <- file.path(output_dir, id)
    run_model(id)
    if (file.exists(out_file)) {
      results_list[[i]] <- mutate(parse_run(out_file), rep = i)
      file.remove(out_file)
    }
    if (file.exists(paste0(out_file, ".netfil"))) file.remove(paste0(out_file, ".netfil"))
  }
} else {
  # --fast: build the population once, then N reps in parallel, each an
  # independent disease-seeding + transmission draw on that shared population.
  # See the WARNING at the top of this file — this mode has been observed to
  # produce false FAILs (d up to 2.4) purely from between-invocation
  # population variance. Only use it for quick iteration, never to conclude
  # anything about whether a change is safe.
  message("--fast mode: results are NOT reliable for a real pass/fail call (see warning at top of script).")
  message("Warm-up run to build population cache...")
  run_model("warmup")
  suppressWarnings(file.remove(file.path(output_dir, "warmup")))
  suppressWarnings(file.remove(file.path(output_dir, "warmup.netfil")))

  plan(multisession, workers = N_WORKERS)
  message(sprintf("Running %d reps in parallel (%d workers)...", N_REPS, N_WORKERS))
  results_list <- future_map(seq_len(N_REPS), function(i) {
    id       <- sprintf("regtest_%04d", i)
    out_file <- file.path(output_dir, id)
    run_model(id)
    if (file.exists(out_file)) {
      out <- mutate(parse_run(out_file), rep = i)
      file.remove(out_file)
      if (file.exists(paste0(out_file, ".netfil"))) file.remove(paste0(out_file, ".netfil"))
      out
    }
  }, .options = furrr_options(seed = TRUE))
  plan(sequential)
}

current <- bind_rows(results_list)
if (nrow(current) == 0) stop("No successful runs — check that model/main runs at all (try running it by hand).")
message(sprintf("Collected %d / %d successful reps.", n_distinct(current$rep), N_REPS))

# ── Update baseline mode: write and stop ─────────────────────────────────────
if (UPDATE_BASELINE) {
  write_csv(current, baseline_csv)
  git_sha <- tryCatch(
    system2("git", c("-C", script_dir, "rev-parse", "--short", "HEAD"), stdout = TRUE),
    error = function(e) "unknown"
  )
  writeLines(c(
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Git commit: %s", git_sha),
    sprintf("Scenario: %s / %s / Theta_%s, no MDA", COUNTRY, SCALE, THETA2_SUFFIX),
    sprintf("T1=%s  K=%s  W=%s  sigma_group=%s  beta_0=%s",
            REF_T1, REF_K, REF_W, REF_SIGMA_GROUP, REF_BETA_0),
    sprintf("n_reps=%d  fast_mode=%s", N_REPS, FAST_MODE)
  ), baseline_meta)
  message(sprintf("\nBaseline written to %s (%d reps). Commit this along with your changes.",
                  baseline_csv, n_distinct(current$rep)))
  quit(status = 0)
}

# ── Compare mode: load baseline, run stats ───────────────────────────────────
if (!file.exists(baseline_csv)) {
  stop("No baseline found at ", baseline_csv, " — run with --update-baseline first.")
}
baseline <- read_csv(baseline_csv, show_col_types = FALSE)

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (pooled_sd == 0) return(0)
  (mean(y) - mean(x)) / pooled_sd
}

verdict_for <- function(d) {
  ad <- abs(d)
  if (ad >= D_FAIL) "FAIL" else if (ad >= D_WARN) "WARN" else "PASS"
}

compare_one <- function(yr, metric_col) {
  b <- baseline |> filter(year == yr, day == 273) |> pull({{ metric_col }})
  c <- current  |> filter(year == yr, day == 273) |> pull({{ metric_col }})
  if (length(b) < 2 || length(c) < 2) return(NULL)
  d  <- cohens_d(b, c)
  tt <- t.test(b, c)
  tibble(
    year = yr, metric = as_label(enquo(metric_col)),
    baseline_mean = mean(b), baseline_sd = sd(b),
    current_mean  = mean(c), current_sd  = sd(c),
    cohens_d = d, p_value = tt$p.value,
    verdict = verdict_for(d)
  )
}

years <- sort(intersect(unique(baseline$year), unique(current$year)))
comparison <- bind_rows(
  map(years, ~ compare_one(.x, ant_prev)),
  map(years, ~ compare_one(.x, mf_prev))
) |> arrange(year, metric)

# ── Report ────────────────────────────────────────────────────────────────
sep <- strrep("=", 70)
cat("\n", sep, "\n", sep = "")
cat(sprintf("  Regression test: current code vs baseline\n"))
cat(sprintf("  Scenario: %s / %s / Theta_%s, no MDA | baseline n=%d, current n=%d\n",
            COUNTRY, SCALE, THETA2_SUFFIX, n_distinct(baseline$rep), n_distinct(current$rep)))
cat(sep, "\n", sep = "")
cat(sprintf("%-6s %-6s %-16s %-16s %10s %10s %8s\n",
            "Year", "Metric", "Baseline", "Current", "Cohen's d", "p", "Verdict"))
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  cat(sprintf("%-6d %-6s %6.3f ± %-6.3f %6.3f ± %-6.3f %10.3f %10.4f %8s\n",
              r$year, r$metric, r$baseline_mean, r$baseline_sd,
              r$current_mean, r$current_sd, r$cohens_d, r$p_value, r$verdict))
}
cat(strrep("-", 70), "\n", sep = "")

overall <- if (any(comparison$verdict == "FAIL")) "FAIL" else if (any(comparison$verdict == "WARN")) "WARN" else "PASS"
worst_row <- comparison |> slice_max(abs(cohens_d), n = 1)
cat(sprintf("OVERALL: %s  (max |d| = %.3f, at year %d %s)\n",
            overall, abs(worst_row$cohens_d[1]), worst_row$year[1], worst_row$metric[1]))
cat(sep, "\n\n", sep = "")

if (overall == "WARN") {
  message("WARN: some effect sizes are in the 'small-to-medium' range. ",
          "Consider re-running with --n=50 or more for a clearer signal before concluding anything.")
}
if (overall == "FAIL") {
  message("FAIL: at least one comparison shows a medium-or-larger effect size. ",
          "This looks like a real behavioral change, not noise — worth investigating ",
          "before merging (see the old-vs-new bisection method used earlier in this ",
          "project's history for how to localize it).")
}

quit(status = if (overall == "FAIL") 1 else 0)
