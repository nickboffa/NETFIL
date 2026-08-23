library(sf)
library(ggplot2)

vil <- st_read("data_storage/wsm_district_boundaries.geojson")
#vil <- st_read("data/tutuila_villages.geojson")

# # below three were merged in 2016 survey paper as a PSU
# merge_set <- c("Satala", "Anua", "Atuu")
# 
# # --- check each target name exists BEFORE merging ---
# setNames(merge_set %in% vil$VILLAGE, merge_set) 
# 
# present <- merge_set[merge_set %in% vil$VILLAGE]
# 
# vil <- vil %>%
#   mutate(psu_name = ifelse(VILLAGE %in% present,
#                            "Satala-Anua-Atuu", VILLAGE)) %>%
#   group_by(psu_name) %>%
#   summarise(.groups = "drop")     # dissolves shared boundaries into one polygon
# 
# st_write(vil, "data/tutuila_psus.geojson")

cent <- st_centroid(st_geometry(vil))
nb   <- st_relate(vil, pattern = "F***1****")   # rook adjacency: shared edge, not just a corner

# mean distance over each unordered neighbour pair
dists <- unlist(lapply(seq_along(nb), function(i) {
  j <- nb[[i]][nb[[i]] > i]                       # avoid double-counting pairs
  if (length(j)) as.numeric(st_distance(cent[i], cent[j])) else numeric(0)
}))

mean(dists)
median(dists)

## PLOTTING

# --- load & project to metres (REQUIRED for distances/labels to be real) ---
vil <- st_transform(vil, 32602)   # UTM zone for American Samoa (EPSG:32602); adjust to your grid's CRS

# --- centroids ---
cent <- st_centroid(st_geometry(vil))
cent_sf <- st_sf(id = seq_len(nrow(vil)), geometry = cent)
xy <- st_coordinates(cent)

# --- rook adjacency: shared edge, not corner-only ---
nb <- st_relate(vil, pattern = "F***1****")

# --- build unique neighbour pairs (i < j to avoid double-counting) ---
pairs <- do.call(rbind, lapply(seq_along(nb), function(i) {
  j <- nb[[i]][nb[[i]] > i]
  if (length(j)) data.frame(i = i, j = j) else NULL
}))

# --- one LINESTRING per pair, with length in metres ---
lines <- st_sfc(
  lapply(seq_len(nrow(pairs)), function(k) {
    st_linestring(rbind(xy[pairs$i[k], ], xy[pairs$j[k], ]))
  }),
  crs = st_crs(vil)
)
lines_sf <- st_sf(pairs, geometry = lines)
lines_sf$dist_m <- round(as.numeric(st_length(lines_sf)))

# --- midpoints for labels ---
mid <- st_coordinates(st_centroid(lines_sf))
lines_sf$mx <- mid[, 1]
lines_sf$my <- mid[, 2]

# --- plot ---
ggplot() +
  geom_sf(data = vil, fill = "grey95", colour = "grey70", linewidth = 0.3) +
  geom_sf(data = lines_sf, colour = "steelblue", linewidth = 0.4) +
  geom_sf(data = cent_sf, colour = "firebrick", size = 1.5) +
  coord_sf(datum = st_crs(vil)) +
  theme_void() +
  labs(x = NULL, y = NULL)

