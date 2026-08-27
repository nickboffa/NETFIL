# ICC -> (sigma_group, beta_0): the group-level random-effect SD implied by
# a target intraclass correlation, and the logit intercept that keeps mean
# prevalence at ant_0 regardless of ICC. Shared by R/write_clustering_params.R
# (ICC derived from a fitted distance-decay curve) and
# R/sequential_abc_raster.R (ICC is a free ABC parameter, passed directly).
icc_to_sigma_beta <- function(icc, ant_0) {
  sigma_group <- sqrt(icc / (1 - icc) * pi^2 / 3)

  expected_prev <- function(mu) {
    integrate(
      function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma_group),
      lower = -Inf, upper = Inf
    )$value
  }

  # expected_prev(mu) is monotonic in mu, so any sufficiently negative lower
  # bound brackets a root — but a fixed -10 isn't negative enough once
  # sigma_group gets large (high ICC), causing uniroot() to fail. Scale the
  # bound with sigma_group so it stays valid across a wide ICC range.
  beta_0 <- uniroot(
    f        = function(mu) expected_prev(mu) - ant_0,
    interval = c(-10 - 10 * sigma_group, 0)
  )$root

  list(sigma_group = sigma_group, beta_0 = beta_0)
}
