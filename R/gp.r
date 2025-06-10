library(ggplot2)
library(ggpubr)
library(dplyr)
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

# the GP will be trained over time and space
gp_data <- pred_data |>
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
