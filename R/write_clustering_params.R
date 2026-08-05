library(jsonlite)

x <- c(20, 1488.961)     # avg distance between polygons, 1488.961 from calc_avg_dist_between_villages median(dists)
y <- c(0.69, 0.30)       # ICCs
fit <- lm(y ~ x)

distance_1 <- 108.256

calc_sigma_and_beta <- function(distance, ant_0) {
  icc <- predict(fit, newdata = data.frame(x = distance))
  sigma_group <- sqrt(icc / (1 - icc) * pi^2 / 3)

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

# Raster scale factor k: folder name is Raster(k*110), distance is k * distance_1
country_scales <- list(
  ASM = c(2, 5, 6),   # Raster220, Raster550, Raster660
  WSM = c(5)          # Raster550
)

for (country in names(country_scales)) {
  cfg   <- fromJSON(file.path("data", country, "country.json"))
  ant_0 <- cfg$seeding$init_ant_prev

  for (k in country_scales[[country]]) {
    params <- calc_sigma_and_beta(k * distance_1, ant_0)
    path   <- file.path("data", country, "Scales", paste0("Raster", k * 110), "clustering_params.csv")
    write.csv(params, path, quote = FALSE, row.names = FALSE)
  }
}
