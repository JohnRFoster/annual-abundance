create_events_per_year <- function(df) {
  df |>
    group_by(propertyID, agrp_prp_id, property_area_km2, year) |>
    summarise(
      all_take = sum(take),
      all_events = n(),
      all_units = sum(trap_count)
    ) |>
    mutate(
      all_take_density = all_take / property_area_km2,
      all_events_density = all_events / property_area_km2,
      all_units_density = all_units / property_area_km2,
      avg_take = all_take / all_events,
      avg_take_density = avg_take / property_area_km2
    )
}

create_methods_per_year <- function(df) {
  df |>
    group_by(propertyID, year, method) |>
    summarise(events = n(), units = sum(trap_count), take = sum(take)) |>
    ungroup() |>
    mutate(
      units_per_event = units / events,
      take_per_unit = take / units,
      take_per_event = take / events
    )
}

create_method_feature <- function(df, col) {
  tmp <- create_methods_per_year(df)

  tmp |>
    select(propertyID, year, method, all_of(col)) |>
    mutate(newname = paste0(method, "_", col)) |>
    select(-method) |>
    pivot_wider(
      names_from = newname,
      values_from = all_of(col),
      values_fill = 0
    )
}

create_all_events_per_year <- function(df) {
  events_per_year <- create_events_per_year(df)
  events_by_method <- create_method_feature(df, "events")
  units_by_method <- create_method_feature(df, "units")
  take_by_method <- create_method_feature(df, "take")
  units_per_event_by_method <- create_method_feature(df, "units_per_event")
  take_per_unit_by_method <- create_method_feature(df, "take_per_unit")
  take_per_event_by_method <- create_method_feature(df, "take_per_event")

  left_join(
    events_per_year,
    events_by_method
  ) |>
    left_join(units_by_method) |>
    left_join(take_by_method) |>
    left_join(units_per_event_by_method) |>
    left_join(take_per_unit_by_method) |>
    left_join(take_per_event_by_method)
}

create_pp_data <- function(df, density_df, first_flag = NA) {
  tmp <- df |>
    group_by(
      propertyID,
      property,
      agrp_prp_id,
      year,
      start_dates,
      end_dates,
      st_name,
      cnty_name,
      farm_bill,
      alws_agrprop_id,
      primary_period,
      property_area_km2,
      county_code,
      trap_count
    ) |>
    summarise(
      total_take_in_pp = sum(take),
      n_events_in_pp = n(),
      n_units_in_pp = sum(trap_count)
    ) |>
    ungroup() |>
    mutate(
      take_density_in_pp = total_take_in_pp / property_area_km2,
      n_events_density_in_pp = n_events_in_pp / property_area_km2,
      n_units_density = n_units_in_pp / property_area_km2,
      farm_bill = if_else(is.na(farm_bill), 0, farm_bill)
    ) |>
    left_join(density_df) |>
    rename(density_estimate = `0.5`) |>
    select(
      propertyID,
      end_dates,
      year,
      st_name,
      cnty_name,
      county_code,
      farm_bill,
      property_area_km2,
      total_take_in_pp,
      take_density_in_pp,
      density_estimate,
      n_events_in_pp,
      n_events_density_in_pp
    )

  tmp |>
    group_by(
      propertyID,
      year,
      st_name,
      county_code,
      property_area_km2,
      farm_bill
    ) |>
    summarise(
      med_density = median(density_estimate),
      all_take = sum(total_take_in_pp),
      avg_take_in_pp = mean(total_take_in_pp),
      all_events = sum(n_events_in_pp),
      avg_events_in_pp = mean(n_events_in_pp),
      n_sampled_pp = n()
    ) |>
    arrange(propertyID, year) |>
    group_by(propertyID) |>
    mutate(
      delta_year = c(first_flag, diff(year)),
      avg_take_density_in_pp = avg_take_in_pp / property_area_km2,
      avg_event_density_in_pp = avg_events_in_pp / property_area_km2,
      events_per_pp = all_events / n_sampled_pp
    ) |>
    ungroup()
}

county_sample_sizes <- function(df) {
  tmp1 <- df |>
    select(propertyID, county_code) |>
    distinct() |>
    count(county_code, name = "n_props_in_county")

  tmp2 <- df |>
    select(propertyID, county_code, year, all_take) |>
    distinct() |>
    group_by(county_code, year) |>
    summarise(n_props_in_county_per_year = n(), county_take = sum(all_take))

  df |>
    left_join(tmp1, by = "county_code") |>
    left_join(tmp2, by = c("county_code", "year"))
}


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
      select(propertyID, st_name, county_code, property_area_km2, farm_bill)

    tmpp <- tibble(
      propertyID = prop_vec[i],
      year = min_year:max_year
    ) |>
      left_join(prop_info, by = "propertyID")

    all_timesteps <- bind_rows(all_timesteps, tmpp)
  }

  left_join(all_timesteps, df) |>
    mutate(
      replace(across(contains("take")), is.na(across(contains("take"))), 0),
      replace(across(contains("events")), is.na(across(contains("events"))), 0),
      n_sampled_pp = replace(n_sampled_pp, is.na(n_sampled_pp), 0)
    )
}

get_take_nn <- function(df, data_repo) {
  lat_lon <- read_csv(file.path(data_repo, "allMISlatLon.csv")) |>
    filter(state_flag == 0)

  tmp <- df |>
    select(propertyID, st_name) |>
    distinct() |>
    left_join(lat_lon) |>
    filter(!is.na(lat), !is.na(lon))

  dist_matrix_m <- tmp |>
    rename(longitude = lon, latititude = lat) |>
    select(longitude, latititude) |>
    geosphere::distm()

  dist_matrix_km <- dist_matrix_m / 1000
  colnames(dist_matrix_km) <- tmp$propertyID

  dist_matrix_km[lower.tri(dist_matrix_km)] <- NA

  nearest_neighbor <- dist_matrix_km |>
    as_tibble() |>
    mutate(propertyID = tmp$propertyID) |>
    pivot_longer(
      cols = -propertyID,
      names_to = "neighbor",
      values_to = "distance_2_nn",
      values_drop_na = TRUE
    ) |>
    filter(distance_2_nn != 0) |>
    group_by(propertyID) |>
    filter(distance_2_nn == min(distance_2_nn)) |>
    mutate(n = seq_len(n())) |>
    filter(n == 1) |>
    ungroup() |>
    select(-n)

  prop_years <- df |>
    select(propertyID, year)

  prop_years_with_neighbor <- left_join(prop_years, nearest_neighbor)

  all_take <- df |>
    select(propertyID, year, contains("take"), -county_take) |>
    rename(neighbor = propertyID)

  nn_with_take <- left_join(prop_years_with_neighbor, all_take) |>
    mutate(
      replace(across(contains("take")), is.na(across(contains("take"))), 0)
    ) |>
    mutate(
      idw_all_take = all_take / distance_2_nn,
      idw_avg_take_in_pp = avg_take_in_pp / distance_2_nn,
      idw_avg_take_density_in_pp = avg_take_density_in_pp / distance_2_nn
    ) |>
    rename(
      nn_all_take = all_take,
      nn_avg_take_in_pp = avg_take_in_pp,
      nn_avg_take_density_in_pp = avg_take_density_in_pp
    )

  nn_with_take
}


create_ml_data <- function(df, per_year_df, land_cover_df) {
  df |>
    left_join(per_year_df) |>
    left_join(land_cover_df) |>
    mutate(
      y = med_density,
      year_fac = as.factor(year),
      farm_bill = as.factor(farm_bill),
      propertyID = factor(propertyID),
      neighbor = factor(neighbor),
      st_name = factor(st_name),
      county_code = factor(county_code),
      state_year = factor(paste(st_name, year_fac)),
      county_year = factor(paste(county_code, year_fac)),
      property_year = factor(paste(propertyID, year_fac))
    ) |>
    select(-med_density, -agrp_prp_id)
}

get_single_year_properties <- function(df) {
  df |>
    group_by(propertyID) |>
    count() |>
    ungroup() |>
    filter(n == 1) |>
    pull(propertyID) |>
    unique()
}

get_multi_year_properties <- function(df) {
  df |>
    filter(!propertyID %in% single_year_properties) |>
    pull(propertyID) |>
    unique()
}

create_m1 <- function(df, single_property_vec, multi_property_vec) {
  tmp1 <- df |>
    filter(propertyID %in% multi_property_vec) |>
    group_by(propertyID) |>
    mutate(
      density_m1 = c(NA, y[1:(n() - 1)]),
      avg_take_m1 = c(NA, avg_take[1:(n() - 1)]),
      all_take_m1 = c(NA, all_take[1:(n() - 1)]),
      events_m1 = c(NA, all_events[1:(n() - 1)])
    ) |>
    ungroup()

  tmp2 <- df |>
    filter(propertyID %in% single_property_vec) |>
    mutate(density_m1 = NA)

  bind_rows(tmp1, tmp2) |>
    mutate(rowID = seq_len(n()))
}

inverse_yj <- function(y, lambda) {
  if (lambda != 0) {
    x <- (1 + y * lambda)^(1 / lambda) - 1
  }

  if (lambda == 0) {
    x <- exp(y) - 1
  }

  x
}


my_recipe <- function(data) {
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = data) |>

    # recommended order from recipes package

    # Impute (does not apply)
    # Handle factor levels
    step_novel(all_nominal_predictors()) |>

    # Individual transformations for skewness and other issues
    step_YeoJohnson(all_outcomes()) |>
    step_YeoJohnson(all_numeric_predictors()) |>

    # Discretize (if needed and if you have no other choice; does not apply)
    # Create dummy variables
    step_unknown(all_nominal_predictors()) |>
    step_dummy(all_nominal_predictors()) |>

    # Create interactions
    step_interact(terms = ~ density_m1:delta_year) |>

    step_interact(terms = ~ all_take:all_events) |>
    step_interact(terms = ~ all_take:all_units) |>
    step_interact(terms = ~ all_events:all_units) |>

    step_interact(terms = ~ all_take_density:all_events_density) |>
    step_interact(terms = ~ all_take_density:all_units_density) |>
    step_interact(terms = ~ all_events_density:all_units_density) |>

    step_interact(terms = ~ avg_take_in_pp:avg_events_in_pp) |>
    step_interact(terms = ~ avg_take_density_in_pp:avg_event_density_in_pp) |>

    step_interact(terms = ~ all_take:rural.road.density) |>
    step_interact(terms = ~ all_take:prop.pub.land) |>
    step_interact(terms = ~ all_take:mean.ruggedness) |>
    step_interact(terms = ~ all_take:mean.canopy.density) |>

    step_interact(terms = ~ avg_take_in_pp:rural.road.density) |>
    step_interact(terms = ~ avg_take_in_pp:prop.pub.land) |>
    step_interact(terms = ~ avg_take_in_pp:mean.ruggedness) |>
    step_interact(terms = ~ avg_take_in_pp:mean.canopy.density) |>

    step_interact(terms = ~ all_take_density:rural.road.density) |>
    step_interact(terms = ~ all_take_density:prop.pub.land) |>
    step_interact(terms = ~ all_take_density:mean.ruggedness) |>
    step_interact(terms = ~ all_take_density:mean.canopy.density) |>

    step_interact(terms = ~ avg_take_density_in_pp:rural.road.density) |>
    step_interact(terms = ~ avg_take_density_in_pp:prop.pub.land) |>
    step_interact(terms = ~ avg_take_density_in_pp:mean.ruggedness) |>
    step_interact(terms = ~ avg_take_density_in_pp:mean.canopy.density) |>

    step_interact(terms = ~ avg_take:rural.road.density) |>
    step_interact(terms = ~ avg_take:prop.pub.land) |>
    step_interact(terms = ~ avg_take:mean.ruggedness) |>
    step_interact(terms = ~ avg_take:mean.canopy.density) |>

    step_interact(terms = ~ avg_take_density:rural.road.density) |>
    step_interact(terms = ~ avg_take_density:prop.pub.land) |>
    step_interact(terms = ~ avg_take_density:mean.ruggedness) |>
    step_interact(terms = ~ avg_take_density:mean.canopy.density) |>

    step_interact(terms = ~ rural.road.density:starts_with("Helicopter")) |>
    step_interact(terms = ~ prop.pub.land:starts_with("Helicopter")) |>
    step_interact(terms = ~ mean.ruggedness:starts_with("Helicopter")) |>
    step_interact(terms = ~ mean.canopy.density:starts_with("Helicopter")) |>

    step_interact(terms = ~ rural.road.density:starts_with("fixedWing")) |>
    step_interact(terms = ~ prop.pub.land:starts_with("fixedWing")) |>
    step_interact(terms = ~ mean.ruggedness:starts_with("fixedWing")) |>
    step_interact(terms = ~ mean.canopy.density:starts_with("fixedWing")) |>

    step_interact(terms = ~ rural.road.density:starts_with("Trap")) |>
    step_interact(terms = ~ prop.pub.land:starts_with("Trap")) |>
    step_interact(terms = ~ mean.ruggedness:starts_with("Trap")) |>
    step_interact(terms = ~ mean.canopy.density:starts_with("Trap")) |>

    step_interact(terms = ~ rural.road.density:starts_with("Snare")) |>
    step_interact(terms = ~ prop.pub.land:starts_with("Snare")) |>
    step_interact(terms = ~ mean.ruggedness:starts_with("Snare")) |>
    step_interact(terms = ~ mean.canopy.density:starts_with("Snare")) |>

    step_interact(terms = ~ rural.road.density:starts_with("Firearms")) |>
    step_interact(terms = ~ prop.pub.land:starts_with("Firearms")) |>
    step_interact(terms = ~ mean.ruggedness:starts_with("Firearms")) |>
    step_interact(terms = ~ mean.canopy.density:starts_with("Firearms")) |>

    # Normalization steps (center, scale, range, etc; does not apply)
    # Multivariate transformation (e.g. PCA, spatial sign, etc; does not apply)

    # filter
    step_nzv(all_predictors())

  blueprint
}
