# Main.R  (Urban air quality - Milestone 3)
library(tidyverse)
library(lubridate)

source("01_download_data.R")

dir.create("outputs", showWarnings = FALSE)

# get urban targets
targets <- download_targets()

# Daily only
targets <- targets %>% filter(duration == "P1D")

targets <- targets %>%
  mutate(datetime = as.POSIXct(datetime, tz = "GMT")) %>%
  arrange(site_id, variable, datetime)

# Visualization only draws a 60-day-window. Is it better?
window_days <- 60  

targets_recent <- targets %>%
  filter(duration == "P1D") %>%
  mutate(datetime = as.POSIXct(datetime, tz = "GMT")) %>%
  filter(datetime >= as.POSIXct(Sys.Date() - window_days, tz = "GMT"))

# get met drivers
site_meta <- download_site_meta()

met <- download_met_drivers(site_meta, past_days = 60) 

# visualization - target time series
p_targets <- targets_recent %>%
  ggplot(aes(x = datetime, y = observation, group = site_id)) +
  geom_line(alpha = 0.6) +
  facet_grid(variable ~ site_id, scales = "free_y") +
  labs(title = "Urban air quality targets (last 60 days)",
       subtitle = paste("Generated:", format(Sys.time(), tz="GMT")),
       x = "Datetime (GMT)", y = "Observation")

ggsave(
  sprintf("outputs/milestone3_targets_last%dd_%s.png", window_days, Sys.Date()),
  p_targets, width = 14, height = 10
)

# visualization - drivers time series
p_met <- met %>%
  ggplot(aes(x = date, y = value, group = site_id)) +
  geom_line(alpha = 0.6) +
  facet_grid(variable ~ site_id, scales = "free_y") +
  labs(title = "Urban meteorological drivers (inputs, example)",
       x = "Date (GMT)", y = "Value")

ggsave("outputs/milestone3_drivers_timeseries.png", p_met, width = 14, height = 8)

message("Milestone 3 plots saved to outputs/.")
