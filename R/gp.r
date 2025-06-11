library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(nimble)

config_name <- "default"
config <- config::get(config = config_name)

out_dir <- config$out_dir
fname <- file.path(out_dir, "allAnnualPredictions.rds")
pred_data <- read_rds(fname)

glimpse(pred_data)

# need lat lon data for properties
data_repo <- config$data_repo
lat_lon <- read_csv(file.path(data_repo, "allMISlatLon.csv")) |>
  filter(state_flag == 0)

# test properties
test_properties <- pred_data |>
  pull(propertyID) |>
  unique() |>
  sample(100)

# the GP will be trained over time and space
gp_data <- pred_data |>
  filter(propertyID %in% test_properties) |>
  select(propertyID, year, pred_density) |>
  left_join(lat_lon, by = "propertyID") |>
  filter(!is.na(lat), !is.na(lon))

# start with all years

make_all_prop_years <- function(df) {
  prop_vec <- df |>
    pull(propertyID) |>
    unique()

  all_timesteps <- tibble()
  for (i in seq_along(prop_vec)) {
    tmp <- df |>
      filter(propertyID == prop_vec[i])

    min_year <- min(tmp$year)
    max_year <- max(tmp$year)

    prop_info <- tmp |>
      dplyr::slice(1) |>
      select(propertyID, st_name, lat, lon)

    tmpp <- tibble(
      propertyID = prop_vec[i],
      year = min_year:max_year
    ) |>
      left_join(prop_info, by = "propertyID")

    all_timesteps <- bind_rows(all_timesteps, tmpp)
  }

  all_timesteps
}

all_years <- make_all_prop_years(gp_data)

time_join <- left_join(all_years, gp_data)

year_cols <- as.character(2014:2023)

time_wide <- time_join |>
  select(-state_flag) |>
  arrange(propertyID, year) |>
  pivot_wider(
    names_from = year,
    values_from = pred_density
  ) |>
  select(propertyID, st_name, lat, lon, all_of(year_cols))

x <- matrix(time_join$year, ncol = 1)
y <- time_join$year
year_matrix <- as.matrix(dist(x, method = "maximum"))
dim(year_matrix)

# normalize the year matrix
year_matrix_normalized <- year_matrix / max(year_matrix)

expcov <- nimbleFunction(
  run = function(
    dists = double(2),
    rho = double(0),
    sigma = double(0),
    tau = double(0)
  ) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) {
          result[i, j] <- pow(sigma, 2) + pow(tau, 2)
        } else {
          result[i, j] <- pow(tau, 2) *
            exp(-rho * pow(dists[i, j], 2))
          result[j, i] <- result[i, j]
        }
      }
    }
    return(result)
  }
)

c_expcov <- compileNimble(expcov)

code <- nimbleCode({
  # priors
  alpha ~ dnorm(0, sd = 100)
  rho ~ dunif(0, 5)
  sigma ~ dunif(0, 100) # prior for variance components based on Gelman (2006)
  tau_cov ~ dunif(0, 100)
  tau_obs ~ dunif(0, 100)

  mu[1:N] <- alpha * ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma, tau_cov)
  s[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])

  # likelihood
  for (i in 1:N) {
    log(mu_x[i]) <- s[i]
    y[i] ~ dnorm(mu_x[i], tau = tau_obs)
  }
})

dists <- as.matrix(year_matrix_normalized)

constants <- list(
  N = nrow(time_join),
  ones = rep(1, nrow(time_join)),
  dists = dists
)

data <- list(
  y = time_join$pred_density
)

inits <- list(
  alpha = rnorm(1),
  sigma = 1,
  tau_obs = 1,
  tau_cov = 1,
  rho = runif(1)
)

## setup initial spatially-correlated latent process values
inits$cov <- c_expcov(
  dists,
  inits$rho,
  inits$sigma,
  inits$tau_cov
)
inits$s <- t(chol(inits$cov)) %*% rnorm(constants$N)
inits$s <- inits$s[, 1] # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
model$initializeInfo()

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
