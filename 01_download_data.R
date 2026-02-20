##' Download Targets
##' @return data.frame in long format with days as rows, and time, site_id, variable, and observed as columns
download_targets <- function(start_dt, end_dt){
  url <- Sys.getenv(
    "URBAN_TARGETS_URL",
    "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/targets/project_id=bu4cast/urban-targets.csv"
  )
  targets <- readr::read_csv(url, col_types = cols())
  
  # Daily targets only - only includes variables PM2.5 and PM10
  targets_recent <- targets %>%
    subset(select = -c(project_id)) %>%
    filter(duration == "P1D") %>%
    mutate(datetime = as.POSIXct(datetime, tz = "GMT")) %>%
    filter(datetime >= start_dt, datetime <= end_dt)
  
  return(targets_recent)
}

##' Download Site metadata
##' @return metadata dataframe
download_site_meta <- function(){
  url <- Sys.getenv(
    "URBAN_SITES_URL",
    "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/targets/project_id=bu4cast/urban-targets-sites.csv"
  )
  readr::read_csv(url, col_types = cols())
}

##' Download met drivers
##' @return metadata dataframe
download_met_drivers <- function(site_meta, past_days) {
  
  # get lat/long
  lat_col <- intersect(names(site_meta), "site_lat")
  lon_col <- intersect(names(site_meta), "site_long")

  lat_col <- lat_col[1]; lon_col <- lon_col[1]
  
  # take a few sites to do milestone 3
  sites_keep <- unique(site_meta$site_id)[1:min(10, n_distinct(site_meta$site_id))]
  site_meta <- site_meta %>% filter(site_id %in% sites_keep)
  
  out <- map_dfr(seq_len(nrow(site_meta)), function(i) {
    site <- site_meta$site_id[i]
    lat  <- site_meta[[lat_col]][i]
    lon  <- site_meta[[lon_col]][i]
    
    resp <- httr::GET(
      "https://api.open-meteo.com/v1/forecast",
      query = list(
        latitude = lat,
        longitude = lon,
        daily = "temperature_2m_mean,precipitation_sum,windspeed_10m_mean",
        past_days = past_days,
        forecast_days = 0,
        timezone = "GMT"
      )
    )
    
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    js  <- jsonlite::fromJSON(txt)
    
    tibble(
      site_id = site,
      date = as.Date(js$daily$time),
      temperature_2m_mean = js$daily$temperature_2m_mean,
      precipitation_sum = js$daily$precipitation_sum,
      windspeed_10m_mean = js$daily$windspeed_10m_mean
    ) %>%
      pivot_longer(-c(site_id, date), names_to = "variable", values_to = "value")
  })
  
  out
}

##' Visualize target time series data
##' @return visualizations in .png
viz_target_time_series <- function(targets, window_days, start_dt, end_dt) {
  p_targets <- targets %>%
    ggplot(aes(x = datetime, y = observation, group = site_id)) +
    geom_line(alpha = 0.6) +
    facet_grid(variable ~ site_id, scales = "free_y") +
    labs(title = sprintf("Urban air quality targets (last %d days)", window_days),
         subtitle = paste0("Data through: ", as.Date(end_dt), " | Generated: ", format(Sys.time(), tz="GMT")),
         x = "Datetime (GMT)", y = "Observation")
  
  ggsave(
    sprintf("outputs/milestone3_targets_last%dd_%s.png", window_days, Sys.Date()),
    p_targets, width = 14, height = 10
  )
}

##' Visualize drivers time series data
##' @return visualizations in .png
viz_drivers_time_series <- function(targets) {
  p_met <- met %>%
    ggplot(aes(x = date, y = value, group = site_id)) +
    geom_line(alpha = 0.6) +
    facet_grid(variable ~ site_id, scales = "free_y") +
    labs(title = "Urban meteorological drivers (inputs, example)",
         x = "Date (GMT)", y = "Value") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ggsave("outputs/milestone3_drivers_timeseries.png", p_met, width = 14, height = 8)
}

# ##' append historical meteorological data into target file
# ##' @param target targets dataframe
# ##' @return updated targets dataframe with added weather data
# merge_met_past <- function(target){
#   
#   ## connect to data
#   df_past <- neon4cast::noaa_stage3()
#   
#   ## filter for site and variable
#   sites <- unique(target$site_id)
#   
#   ## temporary hack to remove a site that's mid-behaving
#   sites = sites[!(sites=="POSE")] 
#   target = target |> filter(site_id %in% sites)  
#   
#   ## grab air temperature from the historical forecast
#   noaa_past <- df_past |> 
#     dplyr::filter(site_id %in% sites,
#                   variable == "air_temperature") |> 
#     dplyr::collect()
#   
#   ## aggregate to daily
#   noaa_past_mean = noaa_past |> 
#     mutate(datetime = as.Date(datetime)) |>
#     group_by(datetime, site_id) |> 
#     summarise(air_temperature = mean(prediction),.groups = "drop")
#   
#   ## Aggregate (to day) and convert units of drivers
#   target <- target %>% 
#     group_by(datetime, site_id,variable) %>%
#     summarize(obs2 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
#     mutate(obs3 = ifelse(is.nan(obs2),NA,obs2)) %>%
#     select(datetime, site_id, variable, obs3) %>%
#     rename(observation = obs3) %>%
#     filter(variable %in% c("temperature", "oxygen")) %>% 
#     tidyr::pivot_wider(names_from = "variable", values_from = "observation")
#   
#   ## Merge in past NOAA data into the targets file, matching by date.
#   target <- left_join(target, noaa_past_mean, by = c("datetime","site_id"))
#   
# }
# 
# ##' Download NOAA GEFS weather forecast
# ##' @param forecast_date start date of forecast
# ##' @return dataframe
# download_met_forecast <- function(forecast_date){
#   noaa_date <- forecast_date - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet
#   
#   ## connect to data
#   df_future <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))
#   
#   ## filter available forecasts by date and variable
#   met_future <- df_future |> 
#     dplyr::filter(datetime >= lubridate::as_datetime(forecast_date), 
#                   variable == "air_temperature") |> 
#     dplyr::collect()
#   
#   ## aggregate to daily
#   met_future <- met_future %>% 
#     mutate(datetime = lubridate::as_date(datetime)) %>% 
#     group_by(datetime, site_id, parameter) |> 
#     summarize(air_temperature = mean(prediction), .groups = "drop") |> 
#     #    mutate(air_temperature = air_temperature - 273.15) |> 
#     select(datetime, site_id, air_temperature, parameter)
#   
#   return(met_future)
# }
