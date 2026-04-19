# 03_fit_and_forecast_pm25_14d.R
# End-to-end: fetch data -> fit state-space model -> use last state as IC -> 14-day ensemble forecast
# Output: forecast dataframe + saved csv/rds in outputs/

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(httr)
  library(jsonlite)
  library(nimble)
  library(coda)
  library(arrow)
})

# config
HIST_DAYS      <- as.integer(Sys.getenv("HIST_DAYS", "365"))       # how many days of history to fit
FORECAST_DAYS  <- as.integer(Sys.getenv("FORECAST_DAYS", "14"))    # horizon
N_ENSEMBLE     <- as.integer(Sys.getenv("N_ENSEMBLE", "100"))      # ensemble members to output
SITE_IDS_ENV   <- Sys.getenv("SITE_IDS", "")                      # optional: comma-separated site_ids
MODEL_ID       <- Sys.getenv("MODEL_ID", "urban_air_forecast")     # used in output
REFERENCE_DATE <- as.Date(Sys.getenv("REFERENCE_DATE", as.character(Sys.Date()))) # reference date
OUT_DIR        <- "outputs/milestone6_like"

dir.create("outputs", showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# source repo functions
if (file.exists("01_download_data.R")) {
  source("01_download_data.R")
} else {
  stop("Cannot find 01_download_data.R in repo root.")
}

# helpper: met download
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

download_met_forecast_one_site <- function(site_meta_one, forecast_days) {
  site <- site_meta_one$site_id[1]

  noaa_forecast_endpoint <- "https://minio-s3.apps.shift.nerc.mghpcc.org"
  noaa_forecast_path <- "bu4cast-ci-read/challenges/project_id=bu4cast/drivers/stage1"
  
  # Today's forecast is probably not ready
  # Yesterday's forecast might not be ready
  # Go back to ereyesterday to be safe
  # (ereyesterday is an actual word that means the day before yesterday)
  today_date <- Sys.Date() |> with_tz("UTC") |> floor_date(unit = "days")
  anchor_date <- today_date - duration(2, "days")
  # Arrow is not clever about date handling; simplify
  anchor_date <- as.character(anchor_date)
  # extend forecast horizon by 2 because our anchor date is ereyesterday
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
    # Temperature, precipitation, and horizontan / vertical wind
    filter(variable %in% c("TMP", "APCP", "VGRD", "UGRD")) |>
    filter(horizon <= forecast_horizon) |>
    # The rest of the filtering can't be done by Arrow :-(
    collect()
  
  if (nrow(stage1_forecast) == 0) {
    stop(paste0("No forecast available for site ", site, " starting on ", anchor_date))
  }
  
  stage1_forecast_cleaned <- stage1_forecast |>
    filter(datetime >= today_date) |>
    pivot_wider(names_from = "variable", values_from = "prediction") |>
    # One forecast per site per day
    mutate(date = floor_date(datetime, unit = "days")) |>
    group_by(site_id, date) |>
    summarize(
      temperature_2m_mean = mean(TMP, na.rm = TRUE),
      precipitation_sum = sum(APCP, na.rm = TRUE),
      windspeed_10m_mean = mean(sqrt(UGRD^2 + VGRD^2), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(stage1_forecast_cleaned)
}

# helper: pick sites
get_sites_to_run <- function(pm25_all) {
  if (nzchar(SITE_IDS_ENV)) {
    str_split(SITE_IDS_ENV, ",")[[1]] %>% str_trim() %>% discard(~ .x == "")
  } else {
    # default: run the top 3 sites with most PM2.5 obs (change to 1 if you want)
    pm25_all %>%
      count(site_id, sort = TRUE, name = "n_obs") %>%
      slice_head(n = 3) %>%
      pull(site_id)
  }
}

# model
build_model_code <- function() {
  nimbleCode({
    # priors (means/sds passed as constants for easy reuse)
    alpha_0     ~ dnorm(pr_mu_alpha0, sd = pr_sd_alpha0)
    alpha_1     ~ dnorm(pr_mu_alpha1, sd = pr_sd_alpha1)
    beta_temp   ~ dnorm(pr_mu_btemp,  sd = pr_sd_btemp)
    beta_wind   ~ dnorm(pr_mu_bwind,  sd = pr_sd_bwind)
    beta_precip ~ dnorm(pr_mu_bprec,  sd = pr_sd_bprec)
    
    # SD priors (kept simple)
    sigma_proc ~ dgamma(2, 2)
    sigma_obs  ~ dgamma(2, 2)
    
    x[1] ~ dnorm(y0, sd = 1)
    y[1] ~ dnorm(x[1], sd = sigma_obs)
    
    for (t in 2:N) {
      mu_x[t] <- alpha_0 +
        alpha_1 * x[t - 1] +
        beta_temp * temp_z[t] +
        beta_wind * wind_z[t] +
        beta_precip * precip_z[t]
      
      x[t] ~ dnorm(mu_x[t], sd = sigma_proc)
      y[t] ~ dnorm(x[t], sd = sigma_obs)
    }
  })
}

default_prior_constants <- function() {
  list(
    pr_mu_alpha0 = 0, pr_sd_alpha0 = 3,
    pr_mu_alpha1 = 0, pr_sd_alpha1 = 1,
    pr_mu_btemp  = 0, pr_sd_btemp  = 1,
    pr_mu_bwind  = 0, pr_sd_bwind  = 1,
    pr_mu_bprec  = 0, pr_sd_bprec  = 1
  )
}

# If you saved param summaries from a previous run, you can use them as informative priors.
# This satisfies "use parameters from previous forecast" in a lightweight way.
load_informative_priors_if_available <- function(site_id) {
  # look for latest param_summary under outputs/
  patt <- paste0("pm25_", site_id, ".*_param_summary\\.csv$")
  files <- list.files("outputs", pattern = patt, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(NULL)
  
  latest <- files[order(file.info(files)$mtime, decreasing = TRUE)][1]
  ps <- readr::read_csv(latest, show_col_types = FALSE)
  
  get_mu_sd <- function(param, default_mu, default_sd) {
    row <- ps %>% filter(parameter == param)
    if (nrow(row) == 0) return(c(default_mu, default_sd))
    mu <- row$median[1]
    # approx sd from 95% interval
    sd <- max((row$q97.5[1] - row$q2.5[1]) / 4, 1e-3)
    c(mu, sd)
  }
  
  pri <- default_prior_constants()
  a0 <- get_mu_sd("alpha_0", pri$pr_mu_alpha0, pri$pr_sd_alpha0)
  a1 <- get_mu_sd("alpha_1", pri$pr_mu_alpha1, pri$pr_sd_alpha1)
  bt <- get_mu_sd("beta_temp", pri$pr_mu_btemp, pri$pr_sd_btemp)
  bw <- get_mu_sd("beta_wind", pri$pr_mu_bwind, pri$pr_sd_bwind)
  bp <- get_mu_sd("beta_precip", pri$pr_mu_bprec, pri$pr_sd_bprec)
  
  list(
    pr_mu_alpha0 = a0[1], pr_sd_alpha0 = a0[2],
    pr_mu_alpha1 = a1[1], pr_sd_alpha1 = a1[2],
    pr_mu_btemp  = bt[1], pr_sd_btemp  = bt[2],
    pr_mu_bwind  = bw[1], pr_sd_bwind  = bw[2],
    pr_mu_bprec  = bp[1], pr_sd_bprec  = bp[2]
  )
}

# ---- fit + forecast for one site ----
fit_and_forecast_one_site <- function(site_id, targets, site_meta) {
  
  site_meta_one <- site_meta %>% filter(site_id == !!site_id)
  if (nrow(site_meta_one) == 0) {
    warning("No site metadata for site_id=", site_id, " ; skipping.")
    return(NULL)
  }
  
  # PM2.5 daily mean
  pm25_site <- targets %>%
    mutate(date = as.Date(datetime)) %>%
    filter(str_detect(variable, regex("PM2\\.?5", ignore_case = TRUE))) %>%
    group_by(site_id, date) %>%
    summarise(pm25 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
    filter(site_id == !!site_id) %>%
    filter(!is.na(pm25)) %>%
    arrange(date)
  
  if (nrow(pm25_site) < 60) {
    warning("Not enough PM2.5 obs for site_id=", site_id, " ; skipping.")
    return(NULL)
  }
  
  # limit history window (last HIST_DAYS)
  end_date <- max(pm25_site$date, na.rm = TRUE)
  start_date <- max(min(pm25_site$date, na.rm = TRUE), end_date - days(HIST_DAYS))
  
  pm25_site <- pm25_site %>% filter(date >= start_date, date <= end_date)
  
  # Met history aligned to the same date window
  met_hist <- download_met_history_one_site(site_meta_one, start_date, end_date)
  
  dat <- pm25_site %>%
    left_join(met_hist, by = c("site_id", "date")) %>%
    arrange(date) %>%
    mutate(
      y = log1p(pm25)
    ) %>%
    drop_na()
  
  if (nrow(dat) < 40) {
    warning("Not enough complete records after joining met for site_id=", site_id, " ; skipping.")
    return(NULL)
  }
  
  # Scale drivers using historical stats
  scale_with <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    if (is.na(sd) || sd == 0) sd <- 1
    list(mu = mu, sd = sd)
  }
  
  s_temp   <- scale_with(dat$temperature_2m_mean)
  s_wind   <- scale_with(dat$windspeed_10m_mean)
  s_precip <- scale_with(log1p(dat$precipitation_sum))
  
  dat <- dat %>%
    mutate(
      temp_z   = (temperature_2m_mean - s_temp$mu) / s_temp$sd,
      wind_z   = (windspeed_10m_mean - s_wind$mu) / s_wind$sd,
      precip_z = (log1p(precipitation_sum) - s_precip$mu) / s_precip$sd
    )
  
  N <- nrow(dat)
  
  # priors: default or from previous run
  pri <- load_informative_priors_if_available(site_id)
  if (is.null(pri)) pri <- default_prior_constants()
  
  code <- build_model_code()
  
  constants <- c(
    list(
      N = N,
      y0 = dat$y[1],
      temp_z = dat$temp_z,
      wind_z = dat$wind_z,
      precip_z = dat$precip_z
    ),
    pri
  )
  
  data_list <- list(y = dat$y)
  
  inits <- function() {
    list(
      alpha_0 = pri$pr_mu_alpha0,
      alpha_1 = pri$pr_mu_alpha1,
      beta_temp = pri$pr_mu_btemp,
      beta_wind = pri$pr_mu_bwind,
      beta_precip = pri$pr_mu_bprec,
      sigma_proc = 1,
      sigma_obs  = 1,
      x = dat$y + rnorm(N, 0, 0.05)
    )
  }
  
  model <- nimbleModel(code = code, constants = constants, data = data_list, inits = inits())
  cmodel <- compileNimble(model)
  
  conf <- configureMCMC(
    model,
    monitors = c("alpha_0","alpha_1","beta_temp","beta_wind","beta_precip","sigma_proc","sigma_obs","x")
  )
  
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  # keep it reasonably fast for daily runs
  samples <- runMCMC(
    cmcmc,
    niter = as.integer(Sys.getenv("MCMC_NITER", "12000")),
    nburnin = as.integer(Sys.getenv("MCMC_NBURN", "3000")),
    thin = as.integer(Sys.getenv("MCMC_THIN", "10")),
    nchains = as.integer(Sys.getenv("MCMC_NCHAINS", "3")),
    samplesAsCodaMCMC = TRUE,
    summary = FALSE,
    WAIC = FALSE
  )
  mat <- do.call(rbind, lapply(samples, as.matrix))
  # pick ensemble draws
  set.seed(1)
  idx <- sample(seq_len(nrow(mat)), size = min(N_ENSEMBLE, nrow(mat)), replace = FALSE)
  draws <- mat[idx, , drop = FALSE]
  
  # last latent state x[N]
  xN_col <- paste0("x[", N, "]")
  if (!xN_col %in% colnames(draws)) stop("Missing last state column: ", xN_col)
  x_curr <- draws[, xN_col]
  
  # drivers forecast (14d) and align dates
  met_fc <- download_met_forecast_one_site(site_meta_one, forecast_days = FORECAST_DAYS + 2) %>%
    arrange(date)
  
  # forecast start = day after last observation date
  fc_start <- end_date + days(1)
  fc_end   <- end_date + days(FORECAST_DAYS)
  
  ### debug ###
  message("end_date = ", end_date)
  message("fc_start = ", fc_start, " | fc_end = ", fc_end)
  message("met_fc nrow = ", nrow(met_fc))
  if (nrow(met_fc) > 0) message("met_fc range: ", min(met_fc$date), " to ", max(met_fc$date))
  ###-------###
  
  # met_fc <- met_fc %>%
  #   filter(date >= fc_start, date <= fc_end)
  # 
  # if (nrow(met_fc) < FORECAST_DAYS) {
  #   warning("Not enough met forecast days returned for site_id=", site_id,
  #           " (got ", nrow(met_fc), ", need ", FORECAST_DAYS, ").")
  # }
  # met_fc already returned forecast starting near today
  met_fc <- met_fc %>%
    arrange(date) %>%
    filter(date >= as.Date(REFERENCE_DATE)) %>%
    slice_head(n = FORECAST_DAYS)
  
  if (nrow(met_fc) < FORECAST_DAYS) {
    warning("Not enough met forecast days for site_id=", site_id,
            " (got ", nrow(met_fc), ", need ", FORECAST_DAYS, ").")
    return(NULL)
  }
  
  # scale forecast drivers using historical scaling
  met_fc <- met_fc %>%
    mutate(
      temp_z   = (temperature_2m_mean - s_temp$mu) / s_temp$sd,
      wind_z   = (windspeed_10m_mean - s_wind$mu) / s_wind$sd,
      precip_z = (log1p(precipitation_sum) - s_precip$mu) / s_precip$sd
    )
  
  H <- nrow(met_fc)
  if (H == 0) return(NULL)
  
  # simulate forward for each ensemble draw
  alpha_0 <- draws[, "alpha_0"]
  alpha_1 <- draws[, "alpha_1"]
  btemp   <- draws[, "beta_temp"]
  bwind   <- draws[, "beta_wind"]
  bprec   <- draws[, "beta_precip"]
  s_proc  <- draws[, "sigma_proc"]
  s_obs   <- draws[, "sigma_obs"]
  
  # matrix: ensemble x horizon
  y_pred <- matrix(NA_real_, nrow = length(idx), ncol = H)
  
  for (e in seq_len(nrow(y_pred))) {
    x_t <- x_curr[e]
    for (h in seq_len(H)) {
      mu_x <- alpha_0[e] +
        alpha_1[e] * x_t +
        btemp[e] * met_fc$temp_z[h] +
        bwind[e] * met_fc$wind_z[h] +
        bprec[e] * met_fc$precip_z[h]
      
      x_t <- rnorm(1, mean = mu_x, sd = s_proc[e])
      y_t <- rnorm(1, mean = x_t, sd = s_obs[e])  # predictive observation on log1p scale
      y_pred[e, h] <- y_t
    }
  }
  
  pm25_pred <- pmax(expm1(y_pred), 0)
  
  # build forecast dataframe (EFI/neon4cast-style)
  reference_dt <- as.POSIXct(REFERENCE_DATE, tz = "GMT")
  
  forecast_df <- expand_grid(
    site_id = site_id,
    datetime = as.POSIXct(met_fc$date, tz = "GMT"),
    ensemble = seq_len(nrow(pm25_pred))
  ) %>%
    mutate(
      reference_datetime = reference_dt,
      variable = "PM2.5",
      family = "ensemble",
      parameter = ensemble,
      prediction = as.vector(t(pm25_pred)),
      model_id = MODEL_ID
    ) %>%
    select(model_id, reference_datetime, datetime, site_id, variable, family, parameter, prediction)
  
  # save artifacts
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "GMT")
  base <- file.path(OUT_DIR, paste0("pm25_", site_id, "_ref_", REFERENCE_DATE, "_", stamp))
  
  saveRDS(samples, paste0(base, "_mcmc.rds"))
  readr::write_csv(dat, paste0(base, "_fit_data.csv"))
  readr::write_csv(forecast_df, paste0(base, "_forecast.csv"))
  
  # lightweight summary (for quick sanity check)
  summ <- forecast_df %>%
    group_by(site_id, datetime) %>%
    summarise(
      p2.5 = quantile(prediction, 0.025, na.rm = TRUE),
      p50  = quantile(prediction, 0.50, na.rm = TRUE),
      p97.5= quantile(prediction, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(summ, paste0(base, "_forecast_summary.csv"))
  
  list(
    site_id = site_id,
    fit_end_date = end_date,
    forecast = forecast_df,
    forecast_summary = summ
  )
}

# ---- main ----
message("Reference date: ", REFERENCE_DATE)
message("History days: ", HIST_DAYS, " | Forecast days: ", FORECAST_DAYS, " | Ensemble: ", N_ENSEMBLE)

# pull latest targets (daily) - give some buffer so we can slice last HIST_DAYS
targets <- download_targets(window_days = max(HIST_DAYS + 30, 120)) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "GMT"))

site_meta <- download_site_meta()

# compute candidate PM2.5 table to decide sites
pm25_all <- targets %>%
  mutate(date = as.Date(datetime)) %>%
  filter(str_detect(variable, regex("PM2\\.?5", ignore_case = TRUE))) %>%
  group_by(site_id, date) %>%
  summarise(pm25 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(pm25))

sites <- get_sites_to_run(pm25_all)
message("Sites to run: ", paste(sites, collapse = ", "))

results <- map(sites, ~ fit_and_forecast_one_site(.x, targets, site_meta))
results <- compact(results)

if (length(results) == 0) stop("No sites produced a forecast (check data availability).")

forecast_all <- bind_rows(map(results, "forecast"))
summ_all <- bind_rows(map(results, "forecast_summary"))

# save combined outputs
combined_base <- file.path(OUT_DIR, paste0("ALL_ref_", REFERENCE_DATE, "_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz="GMT")))
readr::write_csv(forecast_all, paste0(combined_base, "_forecast_ALLSITES.csv"))
readr::write_csv(summ_all, paste0(combined_base, "_forecast_summary_ALLSITES.csv"))

message("Done. Forecast rows: ", nrow(forecast_all))
message("Saved under: ", OUT_DIR)

summ_all |>
  mutate(Site = paste0("Site: ", site_id)) |>
  ggplot(aes(datetime)) +
  geom_ribbon(aes(ymin=`p2.5`, ymax=`p97.5`, fill = "95% CI")) +
  geom_point(aes(y=`p50`, color="PM2.5", shape = "PM2.5"), size = 4) +
  theme_bw(base_size = 16) +
  #coord_cartesian(xlim = c(as.Date("2018-01-01"), as.Date("2019-01-01"))) +
  scale_color_manual(
    name = NULL,
    values = c(
      "PM2.5" = "blue"
    )
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "PM2.5" = "o"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "95% CI" = "lightblue"
    )
  ) +
  ylab(bquote("PM2.5 Concentration" ~ (mu * g %.% m^-3))) +
  xlab("Date") +
  facet_wrap(~Site) +
  ggtitle("Forecast PM2.5 Concentration")
