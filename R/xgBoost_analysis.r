library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)

config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")
source("R/functions_ml.R")

if_dir <- "11_posterior"
top_dir <- config$out_dir

# read in the hyperparameter grid for each eta
etas <- c(0.3, 0.1, 0.05, 0.01, 0.005)
all_ml <- tibble()
for (eta in etas) {
  file_name <- paste0(eta, "ml_YeoJohnsonALL.rds")
  if (!file.exists(file.path(top_dir, file_name))) {
    stop(paste("File not found:", file_name))
  } else {
    out <- read_rds(file.path(top_dir, file_name))
    all_ml <- bind_rows(all_ml, out$hyper_grid)
  }
}

# get the best eta
eta <- all_ml |>
  filter(rmse == min(rmse)) |>
  pull(eta)

best_train <- paste0(eta, "ml_YeoJohnsonALL.rds")

out <- read_rds(file.path(top_dir, best_train))

# fmt: skip
feature_names <- tribble(
  ~Feature,                                ~name,
  "avg_take",                               "Take/event/year",
  "all_take",                               "Take/year",
  "all_events",                             "Events/year",
  "avg_take_x_rural.road.density",          "Take/event/year x county road density",
  "avg_take_x_mean.ruggedness",             "Take/event/year x county ruggedness",
  "avg_take_in_pp_x_mean.ruggedness",       "Take/pp/year x county ruggedness",
  "avg_take_density",                       "Take/event/area/year",
  "avg_take_in_pp",                         "Take/pp/year",
  "avg_take_density_in_pp",                 "Take/area/pp/year",
  "density_m1",                             "Last density estimate",
  "avg_take_density_x_rural.road.density",  "Take/area/year x county road density",
  "avg_take_density_x_mean.ruggedness",     "Take/area/year x county ruggedness",
  "avg_take_density_x_mean.canopy.density", "Take/area/year x county canopy density",
  "density_m1_x_delta_year",                "Last density estimate x years since last estimate",
  "property_area_km2",                      "Property area",
  "county_take",                            "County take",
  "avg_take_in_pp_x_rural.road.density",    "Take/pp/year x county road density",
  "distance_2_nn",                          "Distance (km) to nearest neighbor",
  "mean.ruggedness",                        "County ruggedness",
  "mean.canopy.density",                    "County canopy density"
)

vi_yj <- as_tibble(out$vi) |>
  filter(Feature != "take_per_pp") |>
  dplyr::slice(1:15) |>
  left_join(feature_names)

vi_plot <- function(vi_df) {
  lev <- vi_df |>
    arrange(gainRelative) |>
    pull(name)

  vi_df |>
    mutate(name = factor(name, levels = lev)) |>
    ggplot() +
    aes(x = gainRelative, y = name) +
    geom_point() +
    geom_linerange(aes(xmin = 0, xmax = gainRelative)) +
    labs(
      x = "Relative importance",
      y = "Feature"
    ) +
    theme_bw()
}

vi_plot(vi_yj) +
  labs(title = "Feature importance (YJ transform)") +
  theme(axis.text = element_text(size = 12))

ggsave("plots/ML/variableImportanceAnnualDensity.jpeg")

plot_pred_obs <- function(pred, obs) {
  tibble(
    pred = pred,
    obs = obs
  ) |>
    ggplot() +
    aes(x = obs, y = pred) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 2) +
    # geom_smooth(method = "lm") +
    labs(
      x = "Observed",
      y = "Predicted"
    ) +
    theme_bw()
}

rmsle <- function(pred, obs) {
  sqrt(mean(sum((log(obs + 1) - log(pred + 1))^2)))
}

mean_bias <- function(pred, obs, norm = FALSE) {
  mb <- mean(pred - obs)
  if (norm) {
    return(mb / mean(pred))
  } else {
    return(mb)
  }
}

y <- out$baked_test$y
pred <- out$baked_test$pred
cor(y, pred)
rmsle(y, pred)
mean_bias(y, pred)
mean_bias(y, pred, norm = TRUE)

yj <- plot_pred_obs(pred, y) +
  labs(title = "Predicted-observed (YJ transform)")
yj
ggsave("plots/ML/predObsYJ.jpeg")


lambda <- out$tidy_y$value
y_real <- inverse_yj(y, lambda)
pred_real <- inverse_yj(pred, lambda)

cor(y_real, pred_real)

rmsle(pred_real, y_real)
real <- plot_pred_obs(pred_real, y_real) +
  labs(title = "Predicted-observed (real)")
real
ggsave("plots/ML/predObsReal.jpeg")

ggarrange(yj, real, ncol = 1)
ggsave("plots/ML/predObsBoth.jpeg")

test_names <- tribble(
  ~partition,
  ~name,
  "test_unobserved",
  "Unobserved properties",
  "test_last_year",
  "Last year of properties\nin the training set",
  "test_holdout",
  "New properties"
)

partition_test <- tibble(
  pred = pred,
  obs = y,
  transformation = "YeoJohnson",
  partition = out$baked_test$partition
) |>
  bind_rows(
    tibble(
      pred = pred_real,
      obs = y_real,
      transformation = "None",
      partition = out$baked_test$partition
    )
  ) |>
  left_join(test_names)

partition_test |>
  # filter(obs <= 10) |>
  group_by(partition, transformation) |>
  summarise(
    rmsle = rmsle(pred, obs),
    mbias = mean_bias(pred, obs),
    n_mbias = mean_bias(pred, obs, norm = TRUE),
    cor = cor(pred, obs)
  )

# make two figures

test_plot <- function(df, transform, density_filter = 10000) {
  df |>
    filter(obs <= density_filter, transformation == transform) |>
    ggplot() +
    aes(x = obs, y = pred) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~name, scales = "free_x") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12)
    )
}

g1 <- test_plot(partition_test, "YeoJohnson") +
  labs(
    x = "Observed density (YJ transformed)",
    y = "Predicted density (YJ transformed)"
  )
g2 <- test_plot(partition_test, "None") +
  labs(x = "Observed density (real)", y = "Predicted density (real)") +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggarrange(g1, g2, nrow = 2)

ggsave("plots/ML/predObsTestPartitions.jpeg")

yy <- (((10 + 1)^lambda) - 1) / lambda

g3 <- test_plot(partition_test, "YeoJohnson", yy) +
  labs(
    x = "Observed density (YJ transformed)",
    y = "Predicted density (YJ transformed)"
  )
g4 <- test_plot(partition_test, "None", 10) +
  labs(x = "Observed density (real)", y = "Predicted density (real)") +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggarrange(g3, g4, nrow = 2)
