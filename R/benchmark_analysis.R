# Temporary benchmark analysis — delete after use

library(tidyverse)
library(patchwork)

results <- read_csv("output/benchmark_times.csv") |>
  mutate(
    cp     = as.character(cp),
    pop    = if_else(str_detect(scale, "_\\d+k$"), str_extract(scale, "\\d+k$"), "54k"),
    pop_n  = as.numeric(str_sub(pop, 1, -2)),
    raster = str_remove(scale, "_\\d+k$")
  )

# Summary table
results |>
  group_by(cp, scale, raster, pop, pop_n) |>
  summarise(mean_s = mean(seconds), .groups = "drop") |>
  arrange(cp, raster, pop_n) |>
  print()

# Runtime plot
ggplot(results, aes(x = pop_n, y = seconds, colour = cp, group = cp)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  facet_wrap(~raster) +
  labs(x = "Population (k)", y = "Time (seconds)", colour = "Commuting prop") +
  theme_classic()

# ── Results comparison ────────────────────────────────────────────────────────

output_files <- list.files("output", pattern = "^benchmark_cp.*_rep.*\\.csv$", full.names = TRUE)

epi <- map(output_files, function(f) {
  m <- str_match(basename(f), "benchmark_cp([\\d.]+)_(.+)_rep(\\d+)\\.csv")
  read_csv(f, show_col_types = FALSE) |>
    select(day, year, name, mf_total, ant_total, pop_total) |>
    mutate(
      cp       = m[, 2],
      scale    = m[, 3],
      rep      = m[, 4],
      pop      = if_else(str_detect(scale, "_\\d+k$"), str_extract(scale, "\\d+k$"), "54k"),
      pop_n    = as.numeric(str_sub(pop, 1, -2)),
      raster   = str_remove(scale, "_\\d+k$"),
      mf_prev  = mf_total  / pop_total * 100,
      ant_prev = ant_total / pop_total * 100
    )
}) |>
  list_rbind() |>
  filter(day == 0)

p_mf <- ggplot(epi, aes(x = year, y = mf_prev, colour = as.factor(pop_n),
                         group = interaction(cp, scale, name))) +
  stat_summary(geom = "line", fun = "mean", linewidth = 0.2, color = "black") +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.5, size = 0.2) +
  facet_grid(raster ~ cp) +
  labs(x = "Year", y = "MF prevalence (%)", colour = "Population (k)",
       title = "MF prevalence") +
  theme_classic()

p_ant <- ggplot(epi, aes(x = year, y = ant_prev, colour = as.factor(pop_n),
                          group = interaction(cp, scale, name))) +
  stat_summary(geom = "line", fun = "mean", linewidth = 0.2, color = "black") +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.5, size = 0.2) +
  facet_grid(raster ~ cp) +
  labs(x = "Year", y = "Ag prevalence (%)", colour = "Population (k)",
       title = "Ag prevalence") +
  theme_classic()

(p_mf + p_ant) + plot_layout(guides = "collect")

# Commuting prop

p_mf <- ggplot(epi, aes(x = year, y = mf_prev, colour = cp,
                        group = interaction(cp, raster))) +
  stat_summary(geom = "line", fun = "mean", linewidth = 0.2, color = "black") +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.5, size = 0.2) +
  facet_grid(~raster) +
  labs(x = "Year", y = "MF prevalence (%)", colour = "Cp",
       title = "MF prevalence") +
  theme_classic()

p_ant <- ggplot(epi, aes(x = year, y = ant_prev, colour = as.factor(pop_n),
                         group = interaction(cp, scale, name))) +
  stat_summary(geom = "line", fun = "mean", linewidth = 0.2, color = "black") +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.5, size = 0.2) +
  facet_grid(raster ~ cp) +
  labs(x = "Year", y = "Ag prevalence (%)", colour = "Population (k)",
       title = "Ag prevalence") +
  theme_classic()

(p_mf + p_ant) + plot_layout(guides = "collect")

# ── People per group ──────────────────────────────────────────────────────────

n_groups_per_scale <- tribble(
  ~scale,             ~n_groups,
  "Raster220",        980,
  "Raster220_10k",    957,
  "Raster220_20k",    957,
  "Raster220_100k",   1294,
  "Raster660",        231,
  "Raster660_10k",    230,
  "Raster660_20k",    230,
  "Raster660_100k",   263
)

agents_per_group <- epi |>
  left_join(n_groups_per_scale, by = "scale") |>
  group_by(cp, scale, rep, name) |>
  summarise(mean_agents_per_group = mean(pop_total / n_groups), .groups = "drop")

epi_apg <- epi |>
  left_join(agents_per_group, by = c("cp", "scale", "rep", "name"))

epi_scale <- epi_apg |>
  group_by(cp, scale, year) |>
  summarise(
    mf_prev  = mean(mf_prev),
    ant_prev = mean(ant_prev),
    apg      = mean(mean_agents_per_group),
    .groups  = "drop"
  ) |>
  mutate(cell_size = if_else(str_detect(scale, "660"), 660, 220))

ggplot(epi_scale, aes(year, ant_prev, colour = apg, group = scale)) +
  geom_line(linewidth = 0.6) +
  scale_colour_viridis_c(trans = "log10", name = "Mean agents/group") +
  labs(x = "Year", y = "Ag prevalence (%)", linetype = "Commuting prop") +
  theme_classic()

mean_seconds <- results |>
  group_by(cp, scale) |>
  summarise(seconds = mean(seconds), .groups = "drop")

epi_scale <- epi_scale |>
  left_join(mean_seconds, by = c("scale"))

ggplot(epi_scale, aes(x = apg, y = seconds, colour = cp, group = cp)) +
  geom_point(size = 3) +
  facet_wrap(~cell_size) +
  labs(x = "Mean agents/group", y = "Time (seconds)", colour = "Commuting prop") +
  theme_classic()
