library(tidyverse)

most_recent_forecast <- list.files(
  "outputs/milestone6_like",
  "ALL_ref_.*_forecast_ALLSITES.csv",
  full.names = TRUE
) |>
  sort(decreasing = TRUE) |>
  head(1)

forecast <- read_csv(most_recent_forecast)

uncertainty_time <- forecast %>%
  group_by(datetime) %>%
  summarize(
    mu_t = mean(prediction),
    ensemble_uncertainty = var(prediction),
    .groups = "drop"
  )

var_IC <- var(uncertainty_time$mu_t)
var_ensemble <- mean(uncertainty_time$ensemble_uncertainty)
var_total <- var_IC + var_ensemble

uncertainty_decomp <- tibble(
  component = c("Initial Conditions", "Ensemble Spread"),
  variance = c(var_IC, var_ensemble),
  proportion = c(var_IC, var_ensemble) / var_total
)

print(uncertainty_decomp)

ggplot(uncertainty_time, aes(x = datetime)) +
  geom_line(aes(y = mu_t)) +
  geom_ribbon(
    aes(
      ymin = mu_t - sqrt(ensemble_uncertainty),
      ymax = mu_t + sqrt(ensemble_uncertainty)
    ),
    alpha = 0.3
  ) +
  labs(
    y = "PM2.5",
    title = "Forecast with Uncertainty"
  ) +
  theme_minimal()

V_rel <- uncertainty_time %>%
  mutate(var_total = var_IC + ensemble_uncertainty) %>%
  mutate(
    prop_IC = var_IC / var_total,
    prop_ensemble = ensemble_uncertainty / var_total
  )

ggplot(V_rel, aes(x = datetime)) +
  geom_area(aes(y = prop_IC + prop_ensemble, fill = "Process"), alpha = 0.5) +
  geom_area(aes(y = prop_IC, fill = "IC"), alpha = 0.7) +
  labs(
    title = "Relative Variance Partitioning Over Time",
    y = "Proportion of Variance", x = "Datetime"
  ) +
  scale_fill_manual(
    name = "Source",
    values = list(
      "Process" = "steelblue",
      "IC" = "tomato"
    )
  )
  theme_minimal()
