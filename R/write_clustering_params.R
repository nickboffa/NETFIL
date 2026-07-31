
x <- c(20, 1488.961)     # avg distance between polygons, 1488.961 from calc_avg_dist_between_villages median(dists)
y <- c(0.69, 0.30)       # ICCs
fit <- lm(y ~ x)

distance_1 <- 108.256

ant_0 <- 0.0325

# distance: raster scale
calc_sigma_and_beta <- function(distance) {
  
  icc <- predict(fit, newdata = data.frame(x = c(distance)))
  sigma_group <- sqrt( (icc/(1-icc) * pi^2/3 ) )
  
  ant_0 <- 0.0325
  
  # Expected prevalence under logit-normal(mu, sigma^2)
  expected_prev <- function(mu) {
    integrate(
      function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma_group),
      lower = -Inf, upper = Inf
    )$value
  }
  
  beta_0 <- uniroot(
    f        = function(mu) expected_prev(mu) - ant_0,
    interval = c(-10, 0)
  )$root
  
  data.frame(sigma_group = sigma_group, beta_0 = beta_0)
}

path_base <- "data/Scales/Raster"
for (k in c(2, 5, 6)) {
  distance_k <- k * distance_1
  
  params <- calc_sigma_and_beta(distance_k)
  
  path <- paste0(path_base, k, k, "0/clustering_params.csv")
  write.csv(params, path, quote = F, row.names = F)
}
  


# # Analytic check
# expected_prev(beta_0)  # ≈ 0.0325
# 
# # Monte Carlo cross-check
# mean(plogis(rnorm(1e8, beta_0, sigma_group)))  # ≈ 0.0325
