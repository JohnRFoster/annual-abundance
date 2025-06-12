library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(nimble)

set.seed(920)

source("R/functions_gp.R")
source("R/nimble_gp.R")

config_name <- "default"
config <- config::get(config = config_name)

out_dir <- config$out_dir
fname <- file.path(out_dir, "allAnnualPredictions.rds")
pred_data <- read_rds(fname)

glimpse(pred_data)

# need lat lon data for properties
data_repo <- config$data_repo
lat_lon <- read_csv(file.path(data_repo, "allMISlatLon.csv")) |>
  filter(state_flag == 0) |>
  rename(Lat = lat, Long = lon)

# test properties
test_properties <- pred_data |>
  filter(st_name == "TEXAS") |>
  pull(propertyID) |>
  unique() |>
  sample(50)

# the GP will be trained over time and space
gp_data <- pred_data |>
  filter(propertyID %in% test_properties) |>
  select(propertyID, year, pred_density) |>
  left_join(lat_lon, by = "propertyID") |>
  filter(!is.na(Lat), !is.na(Long))

# start with all years
all_years <- make_all_prop_years_gp(gp_data)

time_join <- left_join(all_years, gp_data)

# create the spatial distance matrix
dist_spatial <- create_spatial_distance_matrix(time_join, normalize = TRUE)

# create the temporal distance matrix
dist_time <- create_year_distance_matrix(time_join, normalize = TRUE)

# check matricies
round(dist_spatial[1:10, 1:10], 2)
round(dist_time[1:10, 1:10], 2)

constants <- list(
  N = nrow(time_join),
  ones = rep(1, nrow(time_join)),
  dist_s = dist_spatial,
  dist_t = dist_time
)

data <- list(
  y = time_join$pred_density
)

inits <- list(
  mu_s = 0,
  mu_t = 0,
  sigma_s = runif(1, 1, 10),
  sigma_t = runif(1, 1, 10),
  tau_s = runif(1),
  tau_t = runif(1),
  rho_s = 0.2,
  rho_t = 0.2,
  tau_obs = runif(1)
)

## setup initial spatially-correlated latent process values
inits$cov_s <- c_expcov(
  dist_spatial,
  inits$rho_s,
  inits$sigma_s,
  inits$tau_s
)
inits$s <- t(chol(inits$cov_s)) %*% rnorm(constants$N)
inits$s <- inits$s[, 1] # so can give nimble a vector rather than one-column matrix

inits$cov_t <- c_expcov(
  dist_time,
  inits$rho_t,
  inits$sigma_t,
  inits$tau_t
)
inits$t <- t(chol(inits$cov_t)) %*% rnorm(constants$N)
inits$t <- inits$t[, 1] # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
model$initializeInfo()
summary(model$mu_x)
str(inits)

c_model <- compileNimble(model)

conf <- configureMCMC(model)
conf$addMonitors("mu_x")

mcmc <- buildMCMC(conf)
c_mcmc <- compileNimble(mcmc)

samples <- runMCMC(
  c_mcmc,
  niter = config$n_iter,
  nburnin = config$n_iter / 2,
  nchains = 3
)

samples_m <- bind_rows(
  as_tibble(as.matrix(samples[[1]])) |> mutate(chain = 1),
  as_tibble(as.matrix(samples[[2]])) |> mutate(chain = 2),
  as_tibble(as.matrix(samples[[3]])) |> mutate(chain = 3)
)

saveRDS(samples_m, file.path(out_dir, "gp_samples.rds"))
saveRDS(time_join, file.path(out_dir, "gp_data.rds"))
saveRDS(constants, file.path(out_dir, "nimble_constants.rds"))
