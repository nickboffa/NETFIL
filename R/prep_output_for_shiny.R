library(tidyverse)
library(sf)

select <- dplyr::select

vil <- read_csv("data/villages_lat_lon.csv")
x_min <<- min(vil$X_Centroid)
y_min <<- min(vil$Y_Centroid)

prep_csv_for_shiny <- function(group_type, filename) {
  out <- read_csv(file.path("output", filename))

  locations <- read_csv(
      file.path("data", "Scales", group_type, "groups.csv")
    ) |>
    mutate(Group = as.character(Group))

  # Compute lat/lon on the small locations table (N groups rows) before the big pivot
  locs_sf <- st_as_sf(
    locations |> mutate(X_raw = X + x_min, Y_raw = Y + y_min),
    coords = c("X_raw", "Y_raw"), crs = 32702
  )
  locs_sf <- st_transform(locs_sf, crs = 4326)
  locations$lon <- st_coordinates(locs_sf)[, "X"]
  locations$lat <- st_coordinates(locs_sf)[, "Y"]

  out |>
    mutate(Time = year + day / 365) |>
    select(sim_i, Time, matches("^mf_\\d+$"), matches("^pop_\\d+$")) |>
    # Pivot mf_ and pop_ together so mf_G is always paired with pop_G
    pivot_longer(
      cols = c(matches("^mf_\\d+$"), matches("^pop_\\d+$")),
      names_to = c(".value", "Group"),
      names_pattern = "([^_]+)_(\\d+)"
    ) |>
    mutate(Mf_prev = 100 * mf / pop) |>
    left_join(locations, "Group")
}
