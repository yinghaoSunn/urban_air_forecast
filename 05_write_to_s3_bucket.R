suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- read in file ---
most_recent_forecast <- list.files(
  "outputs/milestone6_like",
  "ALL_ref_.*_forecast_ALLSITES.csv",
  full.names = TRUE
) |>
  sort(decreasing = TRUE) |>
  head(1)

df <- read_csv(most_recent_forecast)

# ---- determine what to submit ----
# determining whether to submit ensembles or summary stats
hist(df$prediction)
plot(ecdf(df$prediction))
plot(density(df$prediction))
# based on these plots, it seems like we can just submit summary stats and call it a Gamma distribution

# creating summary stats
sum_stats <- df %>%
  group_by(site_id) %>%
  summarise(model_id = first(model_id),
            reference_datetime = first(reference_datetime),
            datetime = first(datetime),
            site_id = first(site_id),
            variable = first(variable),
            mean = mean(prediction),
            sd = sd(prediction),
            conf_interv_02.5 = t.test(prediction, conf.level = 0.975)$conf.int[1],
            conf_interv_97.5 = t.test(prediction, conf.level = 0.975)$conf.int[2],
            .groups = "drop") %>%
  pivot_longer(cols = c(mean, sd, conf_interv_02.5, conf_interv_97.5),
               names_to = "statistic",
               values_to = "value")

# create a temp csv file
temp_file <- tempfile(fileext = ".csv")
write.csv(sum_stats, temp_file, row.names = FALSE)

# ---- submit to S3 bucket ----
# create filename
challenge_name <- 'urban'
team_name <-  'urban_ee585_sp26'
fname <- gsub("outputs/milestone6_like/", "", most_recent_forecast)
url <- paste0("https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-write/challenges/project_id=bu4cast/forecasts/", challenge_name, "/",team_name, "_", fname, "-submission.csv")

# write to bucket
res <- PUT(url, body = upload_file(temp_file))

if (status_code(res) %in% c(200, 204)) {
  message("CSV upload succeeded!")
} else {
  message("Failed. Status: ", status_code(res))
}