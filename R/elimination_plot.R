library(tidyverse)
source("prep_output_for_shiny.R")

all_data <- read_csv("output/r660_2-5MDA_30reps.csv") |> 
  mutate(time = year + day / 365.0)

run_id <- cumsum(c(1, diff(all_data$sim_i) < 0))
all_data <- all_data[run_id == max(run_id), ]

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
  
  postmda_lm <- lm(inf_total ~ time, lm_data)
  time_coef <- postmda_lm$coefficients[2]
  
  end_time <- max(data_i$time)
  end_mf <- data_i[data_i$time == end_time, "inf_total"]
  if (end_mf == 0 || time_coef < 0) {
    TRUE
  } else {
    FALSE
  }
}

eliminated <- sapply(unique(all_data$sim_i), \(i) is_eliminated(all_data, i))
names(eliminated) <- as.character(unique(all_data$sim_i))

elim_df <- data.frame(
  sim_i = unique(all_data$sim_i),
  eliminated = eliminated
)

plot_data <- all_data |>
  group_by(sim_i, n_mda_rounds) |>
  summarise(eliminated = first(eliminated)) |>
  group_by(n_mda_rounds) |>
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

ggplot(plot_data, aes(x=n_mda_rounds, y=prop_eliminated)) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.2) +
  geom_point() +
  geom_line() +
  ggtitle("Raster660") +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  theme_minimal()
