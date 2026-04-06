if (file.exists("03_fit_and_forecast_pm25_14d.R")) {
  source("03_fit_and_forecast_pm25_14d.R")
} else {
  stop("03_fit_and_forecast_pm25_14d.R not found.")
}

if (file.exists("04_Uncertainty_Partitioning.r")) {
  source("04_Uncertainty_Partitioning.r")
} else {
  stop("04_Uncertainty_Partitioning.r not found.")
}

if (file.exists("05_write_to_s3_bucket.R")) {
  source("05_write_to_s3_bucket.R")
} else {
  stop("05_write_to_s3_bucket.R not found.")
}