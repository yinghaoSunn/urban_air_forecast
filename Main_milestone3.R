### Urban air quality forecast workflow ###
# Main.R (Updated for Milestone 3)

## Load required libraries
library(tidyverse)
library(lubridate)

## Load required functions
if (file.exists("01_download_data.R")) source("01_download_data.R")

## Set up directory
dir.create("outputs", showWarnings = FALSE)

## Step 0: Define team name and team members
team_info <- list(
  team_name = "BUSp2026_urban",
  team_list = list(list(
    individualName = list(
      givenName = "Nik",
      surName = "Bates-Haus"
    ),
    organizationName = "Boston University",
    electronicMailAddress = "nikbh@bu.edu"
  ),
  list(
    individualName = list(
      givenName = "Emily",
      surName = "Kim"
    ),
    organizationName = "Boston University",
    electronicMailAddress = "ekim7@bu.edu"
  ),
  list(
    individualName = list(
      givenName = "Radiya",
      surName = "Rafat"
    ),
    organizationName = "Boston University",
    electronicMailAddress = "rrafat@bu.edu"
  ),
  list(
    individualName = list(
      givenName = "Yinghao",
      surName = "Sun"
    ),
    organizationName = "Boston University",
    electronicMailAddress = "sunyh@bu.edu"
  ))
)

## Step 1: Get data
# Take the latest time in the targets data as the end of the window 
# to avoid an empty table due to data not being updated
window_days <- 60 # Only 60 days for easy visualization
end_dt <- max(targets$datetime, na.rm = TRUE)
start_dt <- end_dt - lubridate::days(window_days)

targets <- download_targets(start_dt, end_dt)
site_meta <- download_site_meta()
met <- download_met_drivers(site_meta, past_days = 60)

# visualization - target time series
viz_target_time_series(targets, window_days, start_dt, end_dt)

# visualization - drivers time series
viz_drivers_time_series(targets)

message("Milestone 3 plots saved to outputs/.")
