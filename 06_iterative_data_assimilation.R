# 06_iterative_data_assimilation.R
# Sequential data assimilation 

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(httr)
  library(jsonlite)
  library(arrow)
})

# config
if (file.exists("01_download_data.R")) {
  source("01_download_data.R")
} else if (file.exists("01_download_data (3).R")) {
  source("01_download_data (3).R")
} else {
  stop("Cannot find 01_download_data.R")
}

N_PARTICLES      <- as.integer(Sys.getenv("N_PARTICLES", "200"))
FORECAST_DAYS    <- as.integer(Sys.getenv("FORECAST_DAYS", "14"))
MODEL_ID         <- Sys.getenv("MODEL_ID", "urban_air_forecast")
SITE_IDS_ENV     <- Sys.getenv("SITE_IDS", "")
CALIB_END_DATE   <- as.Date(Sys.getenv("CALIB_END_DATE", "2025-11-30"))
ASSIM_START_DATE <- as.Date(Sys.getenv("ASSIM_START_DATE", "2025-12-01"))
ASSIM_END_DATE   <- as.Date(Sys.getenv("ASSIM_END_DATE", "2025-12-31"))
OUT_DIR          <- Sys.getenv("DA_OUT_DIR", "outputs/data_assimilation")
HIST_DIR         <- Sys.getenv("HIST_DIR", "outputs/milestone5")
SAVE_FORECAST_CSV <- tolower(Sys.getenv("SAVE_FORECAST_CSV", "true")) == "true"
SAVE_ANALYSIS_CSV <- tolower(Sys.getenv("SAVE_ANALYSIS_CSV", "true")) == "true"
##### Debug #####
SKIP_FORECAST <- tolower(Sys.getenv("SKIP_FORECAST", "false")) == "true"
###############

set.seed(as.integer(Sys.getenv("DA_SEED", "1")))

dir.create("outputs", showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# helppers
normalize_log_weights <- function(logw) {
  logw[!is.finite(logw)] <- -1e12
  m <- max(logw)
  w <- exp(logw - m)
  s <- sum(w)
  if (!is.finite(s) || s == 0) return(rep(1 / length(logw), length(logw)))
  w / s
}

systematic_resample <- function(weights) {
  n <- length(weights)
  u0 <- runif(1, 0, 1 / n)
  u  <- u0 + (0:(n - 1)) / n
  cs <- cumsum(weights)
  findInterval(u, cs) + 1L
}

safe_scale <- function(x, mu, sd) {
  if (is.na(sd) || sd == 0) sd <- 1
  (x - mu) / sd
}

summarize_particles <- function(particles, site_id, date, model_id = MODEL_ID) {
  tibble(
    model_id = model_id,
    site_id = site_id,
    datetime = as.POSIXct(date, tz = "GMT"),
    variable = "PM2.5",
    mean = mean(pmax(expm1(particles$x), 0), na.rm = TRUE),
    sd = sd(pmax(expm1(particles$x), 0), na.rm = TRUE),
    q02.5 = quantile(pmax(expm1(particles$x), 0), 0.025, na.rm = TRUE),
    q50   = quantile(pmax(expm1(particles$x), 0), 0.50,  na.rm = TRUE),
    q97.5 = quantile(pmax(expm1(particles$x), 0), 0.975, na.rm = TRUE)
  )
}

particles_to_analysis_df <- function(particles, site_id, reference_date, model_id = MODEL_ID) {
  tibble(
    model_id = model_id,
    reference_datetime = as.POSIXct(reference_date, tz = "GMT"),
    datetime = as.POSIXct(reference_date, tz = "GMT"),
    site_id = site_id,
    variable = "PM2.5",
    family = "ensemble",
    parameter = seq_len(nrow(particles)),
    prediction = pmax(expm1(particles$x), 0)
  )
}

download_met_history_one_site <- function(site_meta_one, start_date, end_date) {
  lat <- site_meta_one$site_lat[1]
  lon <- site_meta_one$site_long[1]
  site <- site_meta_one$site_id[1]

  resp <- httr::GET(
    "https://archive-api.open-meteo.com/v1/archive",
    query = list(
      latitude = lat,
      longitude = lon,
      start_date = as.character(start_date),
      end_date = as.character(end_date),
      daily = "temperature_2m_mean,precipitation_sum,windspeed_10m_mean",
      timezone = "GMT"
    )
  )
  httr::stop_for_status(resp)

  js <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))

  tibble(
    site_id = site,
    date = as.Date(js$daily$time),
    temperature_2m_mean = js$daily$temperature_2m_mean,
    precipitation_sum = js$daily$precipitation_sum,
    windspeed_10m_mean = js$daily$windspeed_10m_mean
  )
}


# does NOT hard-code Sys.Date().
download_met_forecast_one_site <- function(site_meta_one, reference_date, forecast_days) {
  site <- site_meta_one$site_id[1]

  noaa_forecast_endpoint <- "https://minio-s3.apps.shift.nerc.mghpcc.org"
  noaa_forecast_path <- "bu4cast-ci-read/challenges/project_id=bu4cast/drivers/stage1"

  reference_date <- as.Date(reference_date)
  anchor_date <- as.character(reference_date - days(2))
  forecast_horizon <- make_difftime(day = forecast_days + 2)

  stage1_bucket <- arrow::s3_bucket(
    noaa_forecast_path,
    endpoint_override = noaa_forecast_endpoint,
    access_key = Sys.getenv("OSN_KEY"),
    secret_key = Sys.getenv("OSN_SECRET")
  )

  stage1_forecast <- arrow::open_dataset(stage1_bucket) |>
    filter(reference_datetime == anchor_date) |>
    filter(site_id == site) |>
    filter(variable %in% c("TMP", "APCP", "VGRD", "UGRD")) |>
    filter(horizon <= forecast_horizon) |>
    collect()

  if (nrow(stage1_forecast) == 0) {
    stop(paste0(
      "No forecast driver available for site ", site,
      " with anchor date ", anchor_date,
      ". This is the main dependency for the retrospective Dec 1-Dec 31 run."
    ))
  }

  stage1_forecast |>
    mutate(datetime = as.POSIXct(datetime, tz = "GMT")) |>
    filter(as.Date(datetime) >= reference_date) |>
    pivot_wider(names_from = "variable", values_from = "prediction") |>
    mutate(date = as.Date(datetime)) |>
    group_by(site_id, date) |>
    summarise(
      temperature_2m_mean = mean(TMP, na.rm = TRUE),
      precipitation_sum = sum(APCP, na.rm = TRUE),
      windspeed_10m_mean = mean(sqrt(UGRD^2 + VGRD^2), na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(date) |>
    slice_head(n = forecast_days)
}

get_sites_to_run <- function(pm25_all) {
  if (nzchar(SITE_IDS_ENV)) {
    str_split(SITE_IDS_ENV, ",")[[1]] |> str_trim() |> discard(~ .x == "")
  } else {
    pm25_all %>%
      count(site_id, sort = TRUE, name = "n_obs") %>%
      slice_head(n = 3) %>%
      pull(site_id)
  }
}

# posterior helppers
latest_file_or_null <- function(files) {
  if (length(files) == 0) return(NULL)
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

load_latest_historical_mcmc <- function(site_id, hist_dir = HIST_DIR) {
  patt <- paste0("pm25_", site_id, ".*_mcmc\\.rds$")
  f <- latest_file_or_null(list.files(hist_dir, pattern = patt, full.names = TRUE, recursive = TRUE))
  if (is.null(f)) {
    stop("Cannot find historical MCMC for site ", site_id,
         ". Run milestone 5 first or point HIST_DIR to the right folder.")
  }
  readRDS(f)
}

extract_particles_from_historical_mcmc <- function(samples, n_particles, n_obs) {
  mat <- do.call(rbind, lapply(samples, as.matrix))
  keep <- c("alpha_0", "alpha_1", "beta_temp", "beta_wind", "beta_precip", "sigma_proc", "sigma_obs")
  x_col <- paste0("x[", n_obs, "]")
  if (!x_col %in% colnames(mat)) stop("Missing ", x_col, " in historical posterior.")
  need <- c(keep, x_col)
  idx <- sample(seq_len(nrow(mat)), size = min(n_particles, nrow(mat)), replace = FALSE)
  out <- as_tibble(mat[idx, need, drop = FALSE])
  names(out)[names(out) == x_col] <- "x"
  out$weight <- 1 / nrow(out)
  out
}

analysis_state_filename <- function(site_id, reference_date) {
  file.path(OUT_DIR, paste0("analysis_", site_id, "_", as.character(reference_date), ".rds"))
}

analysis_summary_filename <- function(site_id, reference_date) {
  file.path(OUT_DIR, paste0("analysis_", site_id, "_", as.character(reference_date), "_summary.csv"))
}

forecast_filename <- function(site_id, reference_date) {
  file.path(OUT_DIR, paste0("forecast_", site_id, "_ref_", as.character(reference_date), ".csv"))
}

load_previous_analysis_if_available <- function(site_id, reference_date) {
  patt <- paste0("analysis_", site_id, "_(\\d{4}-\\d{2}-\\d{2})\\.rds$")
  files <- list.files(OUT_DIR, pattern = patt, full.names = TRUE)
  if (length(files) == 0) return(NULL)

  dates_chr <- str_match(basename(files), patt)[, 2]
  dates <- as.Date(dates_chr)
  keep <- which(!is.na(dates) & dates < as.Date(reference_date))
  if (length(keep) == 0) return(NULL)

  files <- files[keep]
  dates <- dates[keep]
  f <- files[order(dates, decreasing = TRUE)][1]
  readRDS(f)
}

initialize_particles <- function(site_id, reference_date, dat_site) {
  prev <- load_previous_analysis_if_available(site_id, reference_date)
  if (!is.null(prev)) {
    message("Loaded previous analysis for ", site_id, " before ", reference_date)
    return(prev)
  }

  message("No previous analysis found for ", site_id,
          "; initializing from historical posterior.")
  hist_mcmc <- load_latest_historical_mcmc(site_id)
  extract_particles_from_historical_mcmc(
    samples = hist_mcmc,
    n_particles = N_PARTICLES,
    n_obs = nrow(dat_site)
  )
}


# forecast
update_particles_one_day <- function(particles, obs_row) {
  temp_z   <- obs_row$temp_z[1]
  wind_z   <- obs_row$wind_z[1]
  precip_z <- obs_row$precip_z[1]
  y_obs    <- obs_row$y[1]

  mu_x <- particles$alpha_0 +
    particles$alpha_1 * particles$x +
    particles$beta_temp   * temp_z +
    particles$beta_wind   * wind_z +
    particles$beta_precip * precip_z

  x_pred <- rnorm(nrow(particles), mean = mu_x, sd = particles$sigma_proc)

  if (is.na(y_obs)) {
    particles$x <- x_pred
    particles$weight <- 1 / nrow(particles)
    return(particles)
  }

  logw <- dnorm(y_obs, mean = x_pred, sd = particles$sigma_obs, log = TRUE)
  w <- normalize_log_weights(logw)
  idx <- systematic_resample(w)

  out <- particles[idx, , drop = FALSE]
  out$x <- x_pred[idx]
  out$weight <- 1 / nrow(out)
  out
}

forecast_from_particles <- function(particles, met_fc, site_id, reference_date) {
  H <- nrow(met_fc)
  if (H == 0) return(NULL)

  y_pred <- matrix(NA_real_, nrow = nrow(particles), ncol = H)

  for (e in seq_len(nrow(particles))) {
    x_t <- particles$x[e]
    for (h in seq_len(H)) {
      mu_x <- particles$alpha_0[e] +
        particles$alpha_1[e] * x_t +
        particles$beta_temp[e]   * met_fc$temp_z[h] +
        particles$beta_wind[e]   * met_fc$wind_z[h] +
        particles$beta_precip[e] * met_fc$precip_z[h]

      x_t <- rnorm(1, mean = mu_x, sd = particles$sigma_proc[e])
      y_t <- rnorm(1, mean = x_t, sd = particles$sigma_obs[e])
      y_pred[e, h] <- y_t
    }
  }

  pm25_pred <- pmax(expm1(y_pred), 0)

  expand_grid(
    site_id = site_id,
    datetime = as.POSIXct(met_fc$date, tz = "GMT"),
    ensemble = seq_len(nrow(pm25_pred))
  ) %>%
    mutate(
      reference_datetime = as.POSIXct(reference_date, tz = "GMT"),
      variable = "PM2.5",
      family = "ensemble",
      parameter = ensemble,
      prediction = as.vector(t(pm25_pred)),
      model_id = MODEL_ID
    ) %>%
    select(model_id, reference_datetime, datetime, site_id, variable, family, parameter, prediction)
}

# one-site workflow
run_da_one_site <- function(site_id, targets, site_meta) {
  site_meta_one <- site_meta %>% filter(site_id == !!site_id)
  if (nrow(site_meta_one) == 0) {
    warning("No site metadata for ", site_id)
    return(NULL)
  }

  pm25_site <- targets %>%
    mutate(date = as.Date(datetime)) %>%
    filter(str_detect(variable, regex("PM2\\.?5", ignore_case = TRUE))) %>%
    group_by(site_id, date) %>%
    summarise(pm25 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
    filter(site_id == !!site_id) %>%
    filter(!is.na(pm25)) %>%
    arrange(date)

  if (nrow(pm25_site) < 60) {
    warning("Not enough PM2.5 observations for ", site_id)
    return(NULL)
  }

  hist_start <- min(pm25_site$date)
  hist_end   <- max(pm25_site$date)
  met_hist <- download_met_history_one_site(site_meta_one, hist_start, hist_end)

  dat <- pm25_site %>%
    left_join(met_hist, by = c("site_id", "date")) %>%
    arrange(date) %>%
    mutate(y = log1p(pm25)) %>%
    drop_na()

  if (nrow(dat) < 40) {
    warning("Not enough complete joined records for ", site_id)
    return(NULL)
  }

  # historical scaling used both in update and forecast
  s_temp   <- list(mu = mean(dat$temperature_2m_mean, na.rm = TRUE),
                   sd = sd(dat$temperature_2m_mean, na.rm = TRUE))
  s_wind   <- list(mu = mean(dat$windspeed_10m_mean, na.rm = TRUE),
                   sd = sd(dat$windspeed_10m_mean, na.rm = TRUE))
  s_precip <- list(mu = mean(log1p(dat$precipitation_sum), na.rm = TRUE),
                   sd = sd(log1p(dat$precipitation_sum), na.rm = TRUE))

  dat <- dat %>%
    mutate(
      temp_z   = safe_scale(temperature_2m_mean, s_temp$mu, s_temp$sd),
      wind_z   = safe_scale(windspeed_10m_mean, s_wind$mu, s_wind$sd),
      precip_z = safe_scale(log1p(precipitation_sum), s_precip$mu, s_precip$sd)
    )

  dat_hist <- dat %>% filter(date <= CALIB_END_DATE)
  if (nrow(dat_hist) < 40) {
    stop("Historical window through CALIB_END_DATE is too short for ", site_id)
  }

  particles <- initialize_particles(site_id, ASSIM_START_DATE, dat_hist)

  run_dates <- seq(ASSIM_START_DATE, ASSIM_END_DATE, by = "day")
  out_analysis <- list()
  out_forecast <- list()

  for (ref_date in run_dates) {
    ref_date <- as.Date(ref_date, origin = "1970-01-01")
    message("[", site_id, "] reference_date = ", as.character(ref_date))

    obs_row <- dat %>% filter(date == ref_date)
    if (nrow(obs_row) == 0) {
      warning("No observation row for ", site_id, " on ", ref_date,
              ". Skipping update and carrying particles forward unchanged.")
    } else {
      particles <- update_particles_one_day(particles, obs_row)
    }

    saveRDS(particles, analysis_state_filename(site_id, ref_date))

    analysis_summary <- summarize_particles(particles, site_id, ref_date)
    if (SAVE_ANALYSIS_CSV) {
      write_csv(analysis_summary, analysis_summary_filename(site_id, ref_date))
      write_csv(
        particles_to_analysis_df(particles, site_id, ref_date),
        file.path(OUT_DIR, paste0("analysis_", site_id, "_", ref_date, "_ensemble.csv"))
      )
    }
    out_analysis[[as.character(ref_date)]] <- analysis_summary

    #### Debug #####
    if (!SKIP_FORECAST) {
      met_fc <- download_met_forecast_one_site(
        site_meta_one = site_meta_one,
        reference_date = ref_date,
        forecast_days = FORECAST_DAYS
      )
      
      met_fc <- met_fc %>%
        mutate(
          temp_z   = safe_scale(temperature_2m_mean, s_temp$mu, s_temp$sd),
          wind_z   = safe_scale(windspeed_10m_mean, s_wind$mu, s_wind$sd),
          precip_z = safe_scale(log1p(precipitation_sum), s_precip$mu, s_precip$sd)
        )
      
      fc_df <- forecast_from_particles(
        particles = particles,
        met_fc = met_fc,
        site_id = site_id,
        reference_date = ref_date
      )
      
      if (!is.null(fc_df)) {
        if (SAVE_FORECAST_CSV) write_csv(fc_df, forecast_filename(site_id, ref_date))
        out_forecast[[as.character(ref_date)]] <- fc_df
      }
    }
  }

  list(
    analysis = bind_rows(out_analysis),
    forecast = bind_rows(out_forecast)
  )
}


### main ###
message("CALIB_END_DATE   = ", CALIB_END_DATE)
message("ASSIM_START_DATE = ", ASSIM_START_DATE)
message("ASSIM_END_DATE   = ", ASSIM_END_DATE)
message("FORECAST_DAYS    = ", FORECAST_DAYS)
message("N_PARTICLES      = ", N_PARTICLES)

window_days <- as.integer(max(3650, as.numeric(ASSIM_END_DATE - CALIB_END_DATE) + 400))

targets <- download_targets(window_days = window_days) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "GMT"))

site_meta <- download_site_meta()

pm25_all <- targets %>%
  mutate(date = as.Date(datetime)) %>%
  filter(str_detect(variable, regex("PM2\\.?5", ignore_case = TRUE))) %>%
  group_by(site_id, date) %>%
  summarise(pm25 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(pm25))

sites <- get_sites_to_run(pm25_all)
message("Sites to run: ", paste(sites, collapse = ", "))

results <- map(sites, ~ run_da_one_site(.x, targets, site_meta))
results <- compact(results)

if (length(results) == 0) stop("No sites completed data assimilation.")

analysis_all <- bind_rows(map(results, "analysis"))
forecast_all <- bind_rows(map(results, "forecast"))

write_csv(analysis_all, file.path(OUT_DIR, "ALL_analysis_summary.csv"))

if (nrow(forecast_all) > 0) {
  write_csv(forecast_all, file.path(OUT_DIR, "ALL_forecast_ensemble.csv"))

  # Daily forecast summaries
  forecast_summary_all <- forecast_all %>%
    group_by(model_id, reference_datetime, datetime, site_id, variable) %>%
    summarise(
      mean = mean(prediction, na.rm = TRUE),
      sd = sd(prediction, na.rm = TRUE),
      q02.5 = quantile(prediction, 0.025, na.rm = TRUE),
      q50   = quantile(prediction, 0.50,  na.rm = TRUE),
      q97.5 = quantile(prediction, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(forecast_summary_all, file.path(OUT_DIR, "ALL_forecast_summary.csv"))
} else {
  message("SKIP_FORECAST=true or no forecast rows returned; skipping forecast aggregation.")
}

message("Data assimilation completed.")
message("Outputs written to: ", OUT_DIR)
