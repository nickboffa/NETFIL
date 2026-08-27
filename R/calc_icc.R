library(tidyverse)

icc <- tribble(
  ~level,      ~country, ~icc, ~lo,  ~hi,  ~dist,
  "household", "ASM",    0.69, 0.45, 0.86, 20,
  "household", "WSM",    0.63, 0.34, 0.85, 20,
  "village",   "ASM",    0.30, 0.10, 0.61, 1489,
  "psu",       "WSM",    0.12, 0.01, 0.64, 5000,     # <- compute within-PSU mean pairwise dist
  "region",    "WSM",    0.12, 0.01, 0.64, 40000
) |>
  mutate(
    # logit-scale SE handles the [0,1] boundary better than a raw normal approx
    se_z = (qlogis(hi) - qlogis(lo)) / (2 * qnorm(0.975)),
    w    = 1 / se_z^2
  )

fit_df <- filter(icc, !is.na(dist))

m1 <- nls(icc ~ r0 * exp(-dist / phi),
          data = fit_df, weights = w,
          start = list(r0 = 0.7, phi = 2000))

m2 <- nls(icc ~ c + (r0 - c) * exp(-dist / phi),
          data = fit_df, weights = w,
          start = list(r0 = 0.7, c = 0.1, phi = 1500))


AIC(m1, m2)

predict(m2, newdata = data.frame(dist = 5*distance_1))

plot_fit <- function(mod_fit) {
  xs <- seq(1, 2000)
  ys <- predict(mod_fit, newdata = data.frame(dist = xs))
  plot(xs, ys)
}
