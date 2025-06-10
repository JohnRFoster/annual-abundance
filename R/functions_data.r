# Functions for loading and processing data for the analysis

# for missing values, impute with state mean
mean_impute <- function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Deal with observation covariates, a small percentage of which are missing
get_obs_covars <- function(file) {
  require(dplyr)
  require(readr)

  obs_covs <- file |>
    read_csv(show_col_types = FALSE) |>
    mutate(county_code = sprintf("%05d", FIPS)) |>
    dplyr::select(-starts_with("sd"), -NAME, -FIPS)

  obs_covs <- obs_covs |>
    group_by(STATE_NAME) |>
    mutate(
      rural.road.density = mean_impute(rural.road.density),
      prop.pub.land = mean_impute(prop.pub.land),
      mean.ruggedness = mean_impute(mean.ruggedness),
      mean.canopy.density = mean_impute(mean.canopy.density)
    ) |>
    ungroup() |>
    mutate(
      c_road_den = center_scale(rural.road.density),
      c_rugged = center_scale(mean.ruggedness),
      c_canopy = center_scale(mean.canopy.density)
    )

  targets::tar_assert_true(!any(is.na(obs_covs$c_road_den)))
  targets::tar_assert_true(!any(is.na(obs_covs$c_rugged)))
  targets::tar_assert_true(!any(is.na(obs_covs$c_canopy)))

  obs_covs
}

take_filter <- function(df) {
  zero_take_prp <- df |>
    group_by(propertyID) |>
    summarise(sum_take = sum(take)) |>
    ungroup() |>
    filter(sum_take == 0) |>
    pull(propertyID)

  df |> filter(!propertyID %in% zero_take_prp)
}

# resolve duplicate property values - when there are multiple values, take max
resolve_duplicate <- function(insitu_data) {
  property_areas <- insitu_data |>
    distinct(propertyID, property_area_km2) |>
    group_by(propertyID) |>
    summarize(
      n_areas = length(unique(property_area_km2)),
      property_area_max = max(property_area_km2, na.rm = TRUE),
      # properties with all NA areas get -Inf
      # but the following line changes -Inf to NA
      property_area_km2 = ifelse(
        is.infinite(property_area_max),
        NA,
        property_area_km2
      )
    ) |>
    ungroup()

  insitu_data |>
    left_join(property_areas, by = join_by(propertyID, property_area_km2)) |>
    filter(
      !is.na(property_area_km2),
      property_area_km2 >= 1.8,
      effort > 0
    )
}

create_primary_periods <- function(df, interval, data_repo) {
  require(lubridate)

  end_dates <- unique(sort(df$end.date))
  min_date <- min(end_dates)
  max_date <- max(end_dates)

  start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
  end_dates <- c(start_dates[-1] - 1, max_date)

  targets::tar_assert_identical(length(start_dates), length(end_dates))
  targets::tar_assert_true(min(df$start.date) >= min_date)
  targets::tar_assert_true(max(df$start.date) <= max_date)

  timestep_df <- tibble(start_dates, end_dates) |>
    mutate(timestep = seq_len(n()))
  timestep_df$month <- month(timestep_df$end_dates)
  timestep_df$year <- year(timestep_df$end_dates)

  # for each row in the merged data, insert the integer primary period timestep
  df$timestep <- NA
  message("Assign timesteps...")
  pb <- txtProgressBar(max = nrow(df), style = 1)
  for (i in seq_len(nrow(df))) {
    after_start <- which(timestep_df$start_dates <= df$start.date[i]) |> max()
    before_end <- which(timestep_df$end_dates >= df$end.date[i]) |> min()
    if (after_start == before_end) {
      # then the start and end date is contained within a primary period
      df$timestep[i] <- timestep_df$timestep[before_end]
    } # otherwise, timestep[i] will be left as NA and filtered out later
    setTxtProgressBar(pb, i)
  }
  close(pb)

  write_csv(timestep_df, file.path(data_repo, "timestep_df.csv"))

  df |>
    filter(!is.na(timestep)) |>
    left_join(timestep_df) |>
    arrange(propertyID, timestep)
}

create_timestep_df <- function(df) {
  df |>
    select(propertyID, primary_period) |>
    unique() |>
    arrange(propertyID, primary_period) |>
    group_by(propertyID) |>
    mutate(observed_timestep = seq_len(n())) |>
    ungroup() |>
    mutate(primary_period = primary_period - min(primary_period) + 1)
}

get_last_iteration <- function(path, x) {
  posterior_path <- file.path(path, x)
  if (file.exists(posterior_path)) {
    readRDS(posterior_path)
  } else {
    stop("File not found: ", posterior_path)
  }
}

fix_method_names <- function(df) {
  df |>
    mutate(
      method = if_else(method == "FIREARMS", "Sharpshooting", method),
      method = if_else(method == "FIXED WING", "fixedWing", method),
      method = if_else(method == "HELICOPTER", "Helicopter", method),
      method = if_else(method == "SNARE", "Snare", method),
      method = if_else(method == "TRAPS", "Trap", method)
    )
}
