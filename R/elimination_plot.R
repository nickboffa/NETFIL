library(tidyverse)
#source("prep_output_for_shiny.R")

# "output/r660_2-5MDA_30reps.csv"

all_data <- read_csv("output/mox_ivm_5-8mda_30rep.csv") |> 
  mutate(time = year + day / 365.0)

# run_id <- cumsum(c(1, diff(all_data$sim_i) < 0))
# all_data <- all_data[run_id == max(run_id), ]

# LF was considered to be eliminated in a simulation if, 
# at the end of the simulation period, there were either 
# no m.f. positive people or the slope of a linear 
# least-squares regression on the quarterly territory-level 
# m.f. positive individuals count from 1-year post-MDA until 
# the end of the simulation period was less than or equal to zero.

# full_data: raw csv output from model
is_eliminated <- function(full_data, sim_n) {
  data_i <- full_data |> filter(sim_i == sim_n)
  
  last_mda_year <- unique(data_i$mda_start_year + (data_i$n_mda_rounds-1)*data_i$years_between_rounds)
  
  # MDA performed on day 28, but output is only quarterly. So best to go from day 0
  # Should update MDA to be on day 0 anyway? No, I suspect done this way so that
  # you can easily get the 'start of year'
  lm_data <- data_i |> 
    dplyr::filter(year >= last_mda_year + 1)
  
  postmda_lm <- lm(mf_total ~ time, lm_data)
  time_coef <- postmda_lm$coefficients[2]
  
  end_time <- max(data_i$time)
  end_mf <- data_i[data_i$time == end_time, "mf_total"]
  if (end_mf == 0 || time_coef < 0) {
    TRUE
  } else {
    FALSE
  }
}

get_mf_prev_1y <- function(full_data, sim_n) {
  data_i <- full_data |> filter(sim_i == sim_n)
  
  last_mda_year <- unique(data_i$mda_start_year + (data_i$n_mda_rounds-1)*data_i$years_between_rounds)
  
  data_i |> 
    dplyr::filter(time == last_mda_year + 1) |> 
    mutate(mf_prev = 100*mf_total/pop_total) |> 
    pull(mf_prev)
}

eliminated <- sapply(unique(all_data$sim_i), \(i) is_eliminated(all_data, i))
mf_prev_1y <- sapply(unique(all_data$sim_i), \(i) get_mf_prev_1y(all_data, i))

elim_df <- data.frame(
  sim_i = unique(all_data$sim_i),
  eliminated = eliminated,
  mf_prev_1y = mf_prev_1y
)

all_data <- all_data |>
  left_join(elim_df, by="sim_i")

plot_data <- all_data |>
  group_by(sim_i, name, n_mda_rounds) |>
  summarise(eliminated = first(eliminated)) |>
  group_by(name, n_mda_rounds) |>
  summarise(
    prop_eliminated = mean(eliminated),
    n = n(),
    k = sum(eliminated)
  ) |>
  mutate(
    ci = as.data.frame(Hmisc::binconf(k, n, method = "wilson")),
    ci_lo = ci$Lower,
    ci_hi = ci$Upper
  )

ggplot(plot_data, aes(x=n_mda_rounds, y=prop_eliminated, color=name)) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.2) +
  geom_point() +
  geom_line() +
  ggtitle("Raster660") +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  theme_minimal()



ggplot(df, aes(x = mf_prev, fill = group)) +
  stat_ecdf(geom = "area", alpha = 0.6) +
  labs(x = "mf Prevalence (%)", y = "Density") +
  theme_classic()


ggplot(elim_df, aes(x=mf_prev_1y)) +
  stat_ecdf(geom = "area", alpha = 0.6) +
  xlim(0, 0.2)
