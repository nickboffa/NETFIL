
post <- as.data.frame(fit$adj.values) |> setNames(c("T1", "W", "k"))

prior_long <- prior |> pivot_longer(everything(), names_to = "param", values_to = "value") |> mutate(type = "Prior")
post_long  <- post  |> pivot_longer(everything(), names_to = "param", values_to = "value") |> mutate(type = "Posterior")
combined <- bind_rows(prior_long, post_long)
combined$param <- factor(combined$param, levels = c("T1", "W", "k"))

density_plot <- combined |>
  ggplot(aes(x = value, colour = type, fill = type)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~param, scales = "free") +
  scale_colour_manual(values = c(Prior = "black", Posterior = "red")) +
  scale_fill_manual(values = c(Prior = "black", Posterior = "red")) +
  labs(colour = NULL, fill = NULL) +
  theme_minimal()

density_plot

ggsave("/Users/nicholasboffa/Library/CloudStorage/OneDrive-AustralianNationalUniversity/Uni/2026/lf_modelling/NETFIL/data/Fitted/Raster660/Theta_1/densities.png",
       density_plot, width=10, height=3)



scales <- c("One", "Village", "Raster660", "Many", "Raster220")
base   <- here::here("data", "Fitted")

params <- map_dfr(scales, function(s) {
  read_csv(file.path(base, s, "Theta_1", "TranParams.csv"),
           show_col_types = FALSE) |>
    mutate(Scale = s)
}) |>
  mutate(
    Scale    = factor(Scale, levels = scales),
    WorktoNot = if_else(WorktoNot == 0, NA_real_, WorktoNot)
  ) |>
  pivot_longer(c(Theta_1, Agg, WorktoNot), names_to = "Parameter", values_to = "Value") |>
  mutate(Parameter = recode(Parameter,
                            Theta_1   = "θ1 (transmission)",
                            Agg       = "k (aggregation)",
                            WorktoNot = "w (day/night ratio)"
  ))

params <- params |> 
  mutate(
    n_groups = case_match(
      Scale,
      "One" ~ 1,
      "Village" ~ 64,
      "Many" ~ 316,
      "Raster660" ~ 270,
      "Raster220" ~ 1300
    )
  )

ggplot(params, aes(x = n_groups, y = Value, 
                   label = Scale,
                   group = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  ggrepel::geom_label_repel() +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_minimal() +
  labs(x = "# Groups", y = "Posterior median")


# CV diagnostic plot

cv_res <- cv4abc(
  param   = dplyr::select(complete, T1, W, k),
  sumstat = dplyr::select(complete, Antigen_2016, MF_2016),
  nval    = 200,
  tols    = c(0.005, 0.01, 0.05, 1.0),   # compare multiple tolerances
  method  = "loclinear"
)

plot(cv_res, caption = c("T1", "W", "k"))
summary(cv_res)


# How many particles failed entirely?
message(sprintf("  Complete: %d / %d  (%.1f%% missing)",
                nrow(complete), nrow(particles),
                100 * (1 - nrow(complete)/nrow(particles))))

# Summary of raw summary statistics vs observed targets
cat("\n── Raw particle summary stats ──\n")
print(summary(dplyr::select(complete, Antigen_2016, MF_2016)))
cat(sprintf("  Observed targets: Antigen=%.2f  MF=%.2f\n", OBS_ANT_2016, OBS_MF_2016))

# Are observed targets within the simulated range?
cat(sprintf("  Antigen_2016 range: [%.3f, %.3f]  target=%.2f  %s\n",
            min(complete$Antigen_2016, na.rm=TRUE),
            max(complete$Antigen_2016, na.rm=TRUE),
            OBS_ANT_2016,
            ifelse(OBS_ANT_2016 >= min(complete$Antigen_2016, na.rm=TRUE) &
                     OBS_ANT_2016 <= max(complete$Antigen_2016, na.rm=TRUE), "✓ IN RANGE", "✗ OUT OF RANGE")))
cat(sprintf("  MF_2016 range:      [%.3f, %.3f]  target=%.2f  %s\n",
            min(complete$MF_2016, na.rm=TRUE),
            max(complete$MF_2016, na.rm=TRUE),
            OBS_MF_2016,
            ifelse(OBS_MF_2016 >= min(complete$MF_2016, na.rm=TRUE) &
                     OBS_MF_2016 <= max(complete$MF_2016, na.rm=TRUE), "✓ IN RANGE", "✗ OUT OF RANGE")))

# Distribution of L1 distances
ggplot(tibble(l1 = l1_dist), aes(x = l1)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = quantile(l1_dist, ABC_TOL), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = quantile(l1_dist, ABC_TOL), y = Inf,
           label = sprintf("%.0f%% cutoff", ABC_TOL * 100),
           hjust = -0.1, vjust = 2, color = "red") +
  labs(title = "L1 distance distribution", x = "L1 distance to observed", y = "Count") +
  theme_minimal()

ggplot() +
  geom_point(data = complete, aes(Antigen_2016, MF_2016),
             color = "grey70", alpha = 0.4, size = 1) +
  geom_point(data = accepted, aes(Antigen_2016, MF_2016),
             color = "steelblue", alpha = 0.6, size = 1.5) +
  geom_point(aes(x = OBS_ANT_2016, y = OBS_MF_2016),
             color = "red", shape = 8, size = 4, stroke = 1.5) +
  labs(title = "All particles (grey) vs accepted (blue) vs observed (red)",
       x = "Antigen_2016", y = "MF_2016") +
  theme_minimal()

bind_rows(
  complete |> dplyr::select(T1, W, k) |> mutate(group = "All simulated"),
  accepted |> dplyr::select(T1, W, k) |> mutate(group = "Accepted")
) |>
  pivot_longer(c(T1, W, k), names_to = "param", values_to = "value") |>
  ggplot(aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~param, scales = "free") +
  scale_color_manual(values = c("All simulated" = "grey50", "Accepted" = "steelblue")) +
  scale_fill_manual(values  = c("All simulated" = "grey50", "Accepted" = "steelblue")) +
  labs(title = "Parameter distributions: all simulated vs accepted",
       x = "Value", y = "Density", color = NULL, fill = NULL) +
  theme_minimal()


adj <- as.data.frame(fit$adj.values)
colnames(adj) <- c("T1", "W", "k")

bind_rows(
  prior |> mutate(group = "Prior"),
  adj   |> mutate(group = "Posterior")
) |>
  pivot_longer(c(T1, W, k), names_to = "param", values_to = "value") |>
  ggplot(aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~param, scales = "free") +
  scale_color_manual(values = c("Prior" = "grey50", "Posterior" = "firebrick")) +
  scale_fill_manual(values  = c("Prior" = "grey50", "Posterior" = "firebrick")) +
  labs(title = "Prior vs posterior (regression-adjusted)",
       x = "Value", y = "Density", color = NULL, fill = NULL) +
  theme_minimal()


map_dfr(c(0.01, 0.05, 0.10, 0.25), function(test_tol) {
  acc_test <- complete[l1_dist <= quantile(l1_dist, test_tol), ]
  if (nrow(acc_test) < 20) return(NULL)
  fit_test <- abc(
    target  = obs_vec,
    param   = dplyr::select(acc_test, T1, W, k),
    sumstat = dplyr::select(acc_test, Antigen_2016, MF_2016),
    tol = 1.0, method = "loclinear", transf = "none"
  )
  as.data.frame(fit_test$adj.values) |>
    setNames(c("T1", "W", "k")) |>
    mutate(tol = sprintf("%.0f%%", test_tol * 100))
}) |>
  pivot_longer(c(T1, W, k), names_to = "param", values_to = "value") |>
  ggplot(aes(x = value, color = tol, fill = tol)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~param, scales = "free") +
  labs(title = "Posterior sensitivity to tolerance",
       x = "Value", y = "Density", color = "Tolerance", fill = "Tolerance") +
  theme_minimal()
