library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(rsample)
library(recipes)
library(xgboost)

set.seed(123)

cutoff_date <- ymd("2023-12-31")

config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")
source("R/functions_ml.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file) |>
  select(
    rural.road.density,
    prop.pub.land,
    mean.ruggedness,
    mean.canopy.density,
    county_code
  )

if_dir <- "11_posterior"
top_dir <- config$out_iterative
path <- file.path(top_dir, if_dir)

data_raw <- get_last_iteration(path, "modelData.rds")
data <- data_raw |>
  filter(end_dates <= cutoff_date) |>
  mutate(year = year(end_dates)) |>
  fix_method_names()

density <- get_last_iteration(path, "densitySummaries.rds")

## event metrics ----
all_events_per_year <- create_all_events_per_year(data)

## primary period metrics ----
# (a primary period is a group of events)
yearly_summaries_pp <- create_pp_data(data, density)

## want all timesteps, join
all_timesteps <- make_all_prop_years(yearly_summaries_pp)

with_county_ss <- county_sample_sizes(all_timesteps)
nearest_neighbors <- get_take_nn(with_county_ss, data_repo)

data_joined <- left_join(with_county_ss, nearest_neighbors)

## make ready for xgBoost
data_all <- create_ml_data(data_joined, all_events_per_year, data_obs)

# these properties will become testing properties
# set of properties to draw train/test
all_properties <- unique(data_all$propertyID)
single_year_properties <- get_single_year_properties(data_all)
multi_year_properties <- get_multi_year_properties(data_all)

n_single_year <- length(single_year_properties)
n_multi_year <- length(multi_year_properties)
n_all <- length(all_properties)

assertthat::assert_that(
  n_single_year + n_multi_year == n_all
)

data_with_m1 <- create_m1(
  data_all,
  single_year_properties,
  multi_year_properties
)

y_na <- data_with_m1 |>
  filter(is.na(y)) |>
  mutate(partition = "test_missing", density_m1 = NA)

data_ml <- data_with_m1 |>
  filter(!is.na(y))

take_cols <- grep("take", colnames(data_ml), value = TRUE)
event_cols <- grep("event", colnames(data_ml), value = TRUE)
unit_cols <- grep("unit", colnames(data_ml), value = TRUE)
effort_cols <- c(take_cols, event_cols, unit_cols, "n_sampled_pp")

# keep nearest neighbor and county info
effort_cols <- effort_cols[-grep("nn", effort_cols)]
effort_cols <- effort_cols[-grep("idw", effort_cols)]
effort_cols <- effort_cols[-grep("county", effort_cols)]

# convert all single-year property take data to NA
single_year_df <- data_ml |>
  filter(propertyID %in% single_year_properties) |>
  mutate(
    density_m1 = NA,
    replace(
      across(all_of(effort_cols)),
      !is.na(across(all_of(effort_cols))),
      NA
    )
  )

n_test <- 20
test_draws_multi <- sample.int(n_multi_year, n_test)
test_draws_single <- sample.int(n_single_year, n_test)

# training properties with multiple years without their last year
train1 <- data_ml |>
  filter(propertyID %in% multi_year_properties[-test_draws_multi]) |>
  group_by(propertyID) |>
  filter(year < max(year)) |>
  ungroup()

train2 <- single_year_df |>
  filter(propertyID %in% single_year_properties[-test_draws_single])

df_train_bind <- bind_rows(train1, train2) |>
  mutate(partition = "train")

# testing properties with one year
# are acting as new properties that would have no take data
test1 <- single_year_df |>
  filter(propertyID %in% single_year_properties[test_draws_single]) |>
  mutate(partition = "test_unobserved")

# test properties (last year from training set) with multiple years
test2 <- data_ml |>
  filter(propertyID %in% multi_year_properties[-test_draws_multi]) |>
  group_by(propertyID) |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(partition = "test_last_year")

# held out properties because they only have one year
# density_m1 = NA because we want to test as if these properties never had density estimates
test3 <- data_ml |>
  filter(propertyID %in% multi_year_properties[test_draws_multi]) |>
  mutate(partition = "test_holdout", density_m1 = NA)

df_test_bind <- bind_rows(test1, test2, test3)

df_train <- df_train_bind |>
  filter(!is.na(y))

df_test <- bind_rows(df_test_bind)

assertthat::assert_that(all(!df_test$rowID %in% df_train$rowID))

p_train <- round(nrow(df_train) / nrow(data_ml), 2) * 100
p_test <- round(nrow(df_test) / nrow(data_ml), 2) * 100

message("Train test split across data: ", p_train, ":", p_test)

df_blueprint <- df_train |> select(-partition, -year, -rowID)
blueprint <- my_recipe(df_blueprint)
prepare <- prep(blueprint, training = df_blueprint)
baked_train <- bake(prepare, new_data = df_blueprint)
glimpse(baked_train)

x_train <- baked_train |>
  select(-y) |>
  as.matrix()
y_train <- baked_train |> pull(y)
hist(y_train)

# hyperparameter grid
hyper_grid <- expand_grid(
  eta = config$eta,
  max_depth = 3:8,
  min_child_weight = 0.5,
  subsample = 0.5,
  colsample_bytree = 0.5,

  # pseudo-regularization hyperparameter, controls the complexity of a given tree
  # worth exploring as trees become deeper and when a significant difference between
  # train and test CV error. 0 = no regularization
  gamma = c(0, 1, 10, 100, 1000),

  # L2 regularization (ridge penalty) push coefficients near zero
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000),

  # L1 regularization (lasso penalty) push coefficients all the way to zero
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000),

  # storage
  rmse = 0,
  trees = 0
)

if (config_name == "default") {
  hyper_grid <- hyper_grid[1:10, ]
}

objective <- "reg:squarederror"

n_models <- nrow(hyper_grid)
message("Tuning ", n_models, " models...")

pb <- txtProgressBar(min = 1, max = n_models, style = 1)
for (i in seq_len(n_models)) {
  m <- xgb.cv(
    data = x_train,
    label = y_train,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    verbose = 0,
    params = list(
      eta = hyper_grid$eta[i],
      max_depth = hyper_grid$max_depth[i],
      min_child_weight = hyper_grid$min_child_weight[i],
      subsample = hyper_grid$subsample[i],
      colsample_bytree = hyper_grid$colsample_bytree[i],
      gamma = hyper_grid$gamma[i],
      lambda = hyper_grid$lambda[i],
      alpha = hyper_grid$alpha[i]
    )
  )

  hyper_grid$rmse[i] <- min(m$evaluation_log$test_rmse_mean)
  hyper_grid$trees[i] <- m$best_iteration

  setTxtProgressBar(pb, i)
}
close(pb)

tune_grid <- hyper_grid |>
  filter(rmse == min(rmse))

params <- list(
  eta = pull(tune_grid, eta),
  max_depth = pull(tune_grid, max_depth),
  min_child_weight = pull(tune_grid, min_child_weight),
  subsample = pull(tune_grid, subsample),
  colsample_bytree = pull(tune_grid, colsample_bytree),
  gamma = pull(tune_grid, gamma),
  lambda = pull(tune_grid, lambda),
  alpha = pull(tune_grid, alpha)
)

# train final model
fit <- xgboost(
  params = params,
  data = x_train,
  label = y_train,
  nrounds = pull(tune_grid, trees),
  objective = objective,
  verbose = 0
)

make_prediction <- function(model, new_data) {
  if ("y" %in% colnames(new_data)) {
    new_data <- new_data |> select(-y)
  }

  pred <- predict(fit, as.matrix(new_data), interval = "confidence")

  new_data |> mutate(pred = pred)
}

test2 <- df_test |> select(-partition, -year, -rowID)
baked_test <- bake(prepare, new_data = test2)

df_pred <- make_prediction(fit, baked_test)
df_pred$y <- baked_test$y
df_pred$partition <- df_test$partition

rmse <- sqrt(mean((baked_test$y - df_pred$pred)^2))
cc <- cor(baked_test$y, df_pred$pred)
r2 <- cc^2
vi <- xgb.importance(model = fit)
vi$gainRelative <- vi$Gain / max(vi$Gain)

message("===============================")
message("RMSE: ", round(rmse, 2))
message("COR: ", round(cc, 2))
message("R2: ", round(r2, 2))
message("===============================")

plot(df_pred$y, df_pred$pred)
abline(0, 1)

out_list <- list(
  baked_train = baked_train,
  baked_test = df_pred,
  test_test_rmse = rmse,
  test_test_r2 = r2,
  vi = vi,
  hyper_grid = hyper_grid,
  tidy_y = tidy(prepare, number = 2),
  tidy_x = tidy(prepare, number = 3),
  prepare = prepare,
  fit = fit,
  raw_train = df_train,
  data = bind_rows(df_train, df_test),
  y_na = y_na
)

dest <- config$out_delta
filename <- file.path(dest, paste0(config$eta, "ml_YeoJohnsonALL.rds"))
write_rds(out_list, filename)

fname <- file.path(dest, paste0(config$eta, "_xgb.model"))
xgb.save(fit, fname)
