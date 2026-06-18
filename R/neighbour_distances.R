library(tidyverse)
library(spdep)

scales <- c("Many", "Village", "Raster660")

avg_neighbour_dist <- function(scale) {
  groups <- read_csv(
    here::here("data/Scales", scale, "groups.csv"),
    show_col_types = FALSE
  )

  coords <- as.matrix(select(groups, X, Y))
  nb     <- tri2nb(coords)

  pairs <- map_dfr(seq_len(nrow(groups)), function(i) {
    j_vec <- nb[[i]]
    j_vec <- j_vec[j_vec > i]  # each pair counted once
    if (length(j_vec) == 0) return(tibble())
    tibble(
      i    = i,
      j    = j_vec,
      dist = sqrt((groups$X[i] - groups$X[j_vec])^2 +
                  (groups$Y[i] - groups$Y[j_vec])^2)
    )
  })

  tibble(
    scale             = scale,
    n_groups          = nrow(groups),
    n_neighbour_pairs = nrow(pairs),
    mean_dist_m       = mean(pairs$dist),
    median_dist_m     = median(pairs$dist)
  )
}

results <- map_dfr(scales, avg_neighbour_dist)
print(results, width = Inf)
