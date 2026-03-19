library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(nimble)
library(coda)

if (file.exists("01_download_data.R")) source("01_download_data.R")

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/milestone5", recursive = TRUE, showWarnings = FALSE)

# download historical weather for one site
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
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  js  <- jsonlite::fromJSON(txt)
  
  tibble(
    site_id = site,
    date = as.Date(js$daily$time),
    temperature_2m_mean = js$daily$temperature_2m_mean,
    precipitation_sum = js$daily$precipitation_sum,
    windspeed_10m_mean = js$daily$windspeed_10m_mean
  )
}

#  load targets and choose site

targets <- download_targets(window_days = 3650) %>%
  mutate(date = as.Date(datetime))

pm25_all <- targets %>%
  filter(stringr::str_detect(variable, stringr::regex("PM2\\.?5", ignore_case = TRUE))) %>%
  group_by(site_id, date) %>%
  summarise(pm25 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(pm25))

site_id_env <- Sys.getenv("M5_SITE_ID", "")
if (nzchar(site_id_env)) {
  site_id <- site_id_env
} else {
  site_id <- pm25_all %>%
    count(site_id, sort = TRUE, name = "n_obs") %>%
    slice(1) %>%
    pull(site_id)
}

message("Using site: ", site_id)

site_meta <- download_site_meta() %>%
  filter(site_id == !!site_id)

pm25_site <- pm25_all %>%
  filter(site_id == !!site_id) %>%
  arrange(date)

start_date <- min(pm25_site$date, na.rm = TRUE)
end_date   <- max(pm25_site$date, na.rm = TRUE)

met_site <- download_met_history_one_site(site_meta, start_date, end_date)

# 3. join data
dat <- pm25_site %>%
  left_join(met_site, by = c("site_id", "date")) %>%
  arrange(date) %>%
  mutate(
    y = log1p(pm25),
    temp_z   = as.numeric(scale(temperature_2m_mean)),
    wind_z   = as.numeric(scale(windspeed_10m_mean)),
    precip_z = as.numeric(scale(log1p(precipitation_sum)))
  ) %>%
  select(
    site_id, date, pm25, y,
    temp_z, wind_z, precip_z,
    temperature_2m_mean, windspeed_10m_mean, precipitation_sum
  ) %>%
  drop_na()

if (nrow(dat) < 40) {
  stop("Not enough complete daily records after joining PM2.5 and weather.")
}

N <- nrow(dat)


# 4. Bayesian state-space model
# latent state x[t] is log(PM2.5 + 1)
code <- nimbleCode({
  alpha       ~ dnorm(0, sd = 3)
  phi         ~ dunif(-0.99, 0.99)
  beta_temp   ~ dnorm(0, sd = 1)
  beta_wind   ~ dnorm(0, sd = 1)
  beta_precip ~ dnorm(0, sd = 1)
  
  sigma_proc ~ dunif(0, 3)
  sigma_obs  ~ dunif(0, 3)
  
  x[1] ~ dnorm(y0, sd = 1)
  y[1] ~ dnorm(x[1], sd = sigma_obs)
  
  for (t in 2:N) {
    mu_x[t] <- alpha +
      phi * x[t - 1] +
      beta_temp * temp_z[t] +
      beta_wind * wind_z[t] +
      beta_precip * precip_z[t]
    
    x[t] ~ dnorm(mu_x[t], sd = sigma_proc)
    y[t] ~ dnorm(x[t], sd = sigma_obs)
  }
})

constants <- list(
  N = N,
  y0 = dat$y[1],
  temp_z = dat$temp_z,
  wind_z = dat$wind_z,
  precip_z = dat$precip_z
)

data_list <- list(
  y = dat$y
)

inits <- function() {
  list(
    alpha = 0,
    phi = 0.5,
    beta_temp = 0,
    beta_wind = 0,
    beta_precip = 0,
    sigma_proc = 0.2,
    sigma_obs = 0.2,
    x = dat$y + rnorm(N, 0, 0.05)
  )
}

model <- nimbleModel(
  code = code,
  constants = constants,
  data = data_list,
  inits = inits()
)

cmodel <- compileNimble(model)

conf <- configureMCMC(
  model,
  monitors = c(
    "alpha", "phi", "beta_temp", "beta_wind", "beta_precip",
    "sigma_proc", "sigma_obs", "x"
  )
)

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

samples <- runMCMC(
  cmcmc,
  niter = 25000,
  nburnin = 5000,
  thin = 10,
  nchains = 4,
  samplesAsCodaMCMC = TRUE,
  summary = FALSE,
  WAIC = FALSE
)

# 5. save outputs
prefix <- file.path(
  "outputs", "milestone5",
  paste0("pm25_", site_id, "_", Sys.Date())
)

saveRDS(samples, paste0(prefix, "_mcmc.rds"))
write_csv(dat, paste0(prefix, "_fit_data.csv"))

mat <- do.call(rbind, lapply(samples, as.matrix))

monitor_pars <- c(
  "alpha", "phi", "beta_temp", "beta_wind", "beta_precip",
  "sigma_proc", "sigma_obs"
)

param_summary <- tibble(
  parameter = monitor_pars,
  mean   = colMeans(mat[, monitor_pars, drop = FALSE]),
  q2.5   = apply(mat[, monitor_pars, drop = FALSE], 2, quantile, probs = 0.025),
  median = apply(mat[, monitor_pars, drop = FALSE], 2, quantile, probs = 0.5),
  q97.5  = apply(mat[, monitor_pars, drop = FALSE], 2, quantile, probs = 0.975),
  ess    = as.numeric(coda::effectiveSize(samples[, monitor_pars])),
  rhat   = as.numeric(coda::gelman.diag(samples[, monitor_pars], autoburnin = FALSE)$psrf[, 1])
)

write_csv(param_summary, paste0(prefix, "_param_summary.csv"))

x_cols <- grep("^x\\[", colnames(mat), value = TRUE)
x_cols <- x_cols[order(as.integer(stringr::str_extract(x_cols, "\\d+")))]

state_summary <- tibble(
  date = dat$date,
  observed_pm25 = dat$pm25,
  fitted_lo  = pmax(expm1(apply(mat[, x_cols, drop = FALSE], 2, quantile, probs = 0.025)), 0),
  fitted_med = pmax(expm1(apply(mat[, x_cols, drop = FALSE], 2, quantile, probs = 0.5)), 0),
  fitted_hi  = pmax(expm1(apply(mat[, x_cols, drop = FALSE], 2, quantile, probs = 0.975)), 0)
)

write_csv(state_summary, paste0(prefix, "_state_summary.csv"))

message("Milestone 5 model fit completed.")
message("Site used: ", site_id)
message("Files saved under outputs/milestone5/")