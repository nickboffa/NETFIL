library(tidyverse)

tran     <- read_csv("data/Fitted/Many/Theta_1/TranParams.csv")
theta1   <- scan("data/Fitted/Many/Theta_1/Theta1.txt", quiet = TRUE)
agg      <- scan("data/Fitted/Many/Theta_1/Agg.txt",    quiet = TRUE)
work     <- scan("data/Fitted/Many/Theta_1/Work.txt",   quiet = TRUE)

mode_est <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

summary_stats <- tibble(
  stat      = c("mean", "median", "mode", "median_rounded_4sf", "median_rounded_3sf"),
  Theta_1   = c(mean(theta1),        median(theta1),        mode_est(theta1),
                signif(median(theta1), 4), signif(median(theta1), 3)),
  Agg       = c(mean(agg),           median(agg),           mode_est(agg),
                signif(median(agg), 4),    signif(median(agg), 3)),
  WorktoNot = c(mean(work),          median(work),          mode_est(work),
                signif(median(work), 4),   signif(median(work), 3))
)

target <- tibble(
  stat      = "TranParams",
  Theta_1   = tran$Theta_1,
  Agg       = tran$Agg,
  WorktoNot = tran$WorktoNot
)

bind_rows(summary_stats, target) |> print(n = Inf) |> View()






library(tidyverse)

tran <- read_csv("data/Fitted/Many/Theta_1/TranParams.csv")
fit  <- read_tsv("data/Fitted/Many/Theta_1/fit_s1.tsv")

# Posterior draws are unique (T1, W, k) combos — one row per accepted draw
params <- fit |> distinct(T1, k, Ratio_2014)

mode_est <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

summary_stats <- tibble(
  stat      = c("mean", "median", "mode", "median_rounded_4sf", "median_rounded_3sf"),
  Theta_1   = c(mean(params$T1),        median(params$T1),        mode_est(params$T1),
                signif(median(params$T1), 4), signif(median(params$T1), 3)),
  Agg       = c(mean(params$k),         median(params$k),         mode_est(params$k),
                signif(median(params$k), 4),  signif(median(params$k), 3)),
  WorktoNot = c(mean(params$Ratio_2014), median(params$Ratio_2014), mode_est(params$Ratio_2014),
                signif(median(params$Ratio_2014), 4), signif(median(params$Ratio_2014), 3))
)

# Target values from TranParams.csv
target <- tibble(
  stat      = "TranParams",
  Theta_1   = tran$Theta_1,
  Agg       = tran$Agg,
  WorktoNot = tran$WorktoNot
)

bind_rows(summary_stats, target) |> print(n = Inf)

