library(jsonlite)
source("R/calc_icc.R")
source("R/icc_to_sigma_beta.R")

distance_1 <- 108.256

calc_sigma_and_beta <- function(fit, distance, ant_0) {
  icc <- predict(fit, newdata = data.frame(dist = distance)) |> round(2)
  # icc/ant_0 are cosmetic only, same as sequential_abc_raster.R's
  # clustering_params.csv output — model/initialise_pop.cpp's
  # read_parameters() only reads the first two columns (sigma_group,
  # beta_0) and ignores the rest.
  cbind(as.data.frame(icc_to_sigma_beta(icc, ant_0)), icc = icc, ant_0 = ant_0)
}

# Raster scale factor k: folder name is Raster(k*110), distance is k * distance_1
country_scales <- list(
  #ASM = c(2, 5, 6),   # Raster220, Raster550, Raster660
  WSM = c(5)          # Raster550
)

for (country in names(country_scales)) {
  cfg   <- fromJSON(file.path("data", country, "country.json"))
  ant_0 <- cfg$seeding$init_ant_prev

  for (k in country_scales[[country]]) {
    params <- calc_sigma_and_beta(m2, k * distance_1, ant_0)
    path   <- file.path("data", country, "Scales", paste0("Raster", k * 110), "clustering_params.csv")
    write.csv(params, path, quote = FALSE, row.names = FALSE)
  }
}
