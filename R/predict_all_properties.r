#---------
#
# Workflow for predicting yearly abundances across all MIS
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(recipes)
library(lubridate)
library(ggplot2)
library(ggpubr)

cutoff_date <- ymd("2023-12-31")

# config_name <- "hpc_dev"
config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")
source("R/functions_ml.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo
file <- file.path(data_repo, config$file_mis)
all_take <- read_csv(file, show_col_types = FALSE) |>
	filter(start.date >= lubridate::ymd("2014-01-01")) |>
	mutate(
		cnty_name = if_else(
			grepl("ST ", cnty_name),
			gsub("ST ", "ST. ", cnty_name),
			cnty_name
		),
		cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name)
	)

data_farm_bill <- read_csv(file.path(
	data_repo,
	"All_FB_Agreements_long_2024-05-30.csv"
))

farm_bill_properties <- data_farm_bill |>
	rename(alws_agrprop_id = propertyID) |>
	select(-agreement_name, -property_name) |>
	mutate(farm_bill = 1)

data_mis <- all_take |>
	mutate(property_area_km2 = round(property.size * 0.00404686, 2)) |>
	filter(property_area_km2 >= 1.8, st_name != "HAWAII") |>
	mutate(
		effort = if_else(
			cmp_name %in% c("TRAPS, CAGE", "SNARE"),
			cmp.days,
			cmp.hours
		),
		effort_per = effort / cmp.qty,
		cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)
	) |>
	rename(method = cmp_name, trap_count = cmp.qty) |>
	select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
	distinct() |>
	mutate(propertyID = paste0(agrp_prp_id, "-", alws_agrprop_id)) |>
	arrange(propertyID, start.date, end.date)

data_timestep <- create_primary_periods(data_mis, 4, data_repo) |>
	resolve_duplicate() |> # resolve duplicate property areas
	take_filter() |>
	rename(statefp = st_gsa_state_cd, countyfp = cnty_gsa_cnty_cd) |>
	mutate(
		countyfp = sprintf("%03d", countyfp),
		countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
		county_code = as.numeric(paste0(statefp, countyfp)),
		county_code = sprintf("%05d", county_code)
	) |>
	select(
		propertyID,
		agrp_prp_id,
		alws_agrprop_id,
		start_dates,
		end_dates,
		st_name,
		cnty_name,
		county_code,
		method,
		trap_count,
		take,
		property_area_km2,
		effort,
		effort_per,
		timestep
	) |>
	rename(primary_period = timestep)

# now we have two columns for time
# primary_period is how [interval] sequences are aligned across the data set
# timestep is the sequence of primary periods within a property
timestep_df <- create_timestep_df(data_timestep)

data_pp <- left_join(
	data_timestep,
	timestep_df,
	by = join_by(propertyID, primary_period)
) |>
	mutate(
		primary_period = primary_period - min(primary_period) + 1,
		year = year(end_dates)
	) |>
	filter(end_dates <= cutoff_date) |>
	left_join(farm_bill_properties) |>
	mutate(
		farm_bill = if_else(is.na(farm_bill), 0, 1),
		property = as.numeric(as.factor(propertyID)),
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "fixedWing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method)
	)

## event metrics ----
all_events_per_year <- create_all_events_per_year(data_pp)

density_df <- tibble(`0.5` = NA, propertyID = NA, end_dates = NA)
yearly_summaries_pp <- create_pp_data(data_pp, density_df)

# do not want all timesteps,
# test partition on unobserved properties was not good
# all_timesteps <- make_all_prop_years(yearly_summaries_pp)

with_county_ss <- county_sample_sizes(yearly_summaries_pp)
nearest_neighbors <- get_take_nn(with_county_ss, data_repo)

data_joined <- left_join(with_county_ss, nearest_neighbors)

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

## make ready for xgBoost
data_all <- create_ml_data(data_joined, all_events_per_year, data_obs)

if_dir <- "11_posterior"
iterative_dir <- config$out_iterative

posterior_path <- file.path(iterative_dir, if_dir, "modelData.rds")
data_post_raw <- read_rds(posterior_path)

data <- data_post_raw |>
	filter(end_dates <= cutoff_date) |>
	mutate(
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "fixedWing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method),
		year = year(end_dates)
	)

posterior_path <- file.path(iterative_dir, if_dir, "densitySummaries.rds")
density <- read_rds(posterior_path)

bayes_estimates <- data |>
	select(
		propertyID,
		agrp_prp_id,
		property,
		primary_period,
		end_dates,
		year,
		property_area_km2
	) |>
	distinct() |>
	left_join(density) |>
	group_by(propertyID, year) |>
	summarise(y = median(`0.5`)) |>
	ungroup()

data_with_density <- left_join(data_all, bayes_estimates)

single_year_properties <- get_single_year_properties(data_with_density)
multi_year_properties <- get_multi_year_properties(data_with_density)

data_ml <- create_m1(
	data_with_density,
	single_year_properties,
	multi_year_properties
) |>
	mutate(propertyID = factor(propertyID))

length(unique(data_all$propertyID))
length(unique(data_ml$propertyID))

best_train <- paste0(0.005, "ml_YeoJohnsonALL.rds")
out <- read_rds(file.path(config$out_dir, best_train))
prepare <- out$prepare
baked_test <- bake(prepare, new_data = data_ml)
fit <- out$fit

newdata <- baked_test |>
	select(-y) |>
	as.matrix()

pred <- predict(fit, newdata, interval = "confidence")

lambda <- out$tidy_y$value

data_pred <- data_ml |>
	select(-y) |>
	mutate(pred_yj = pred, pred_density = inverse_yj(pred, lambda))

fname <- file.path(config$out_dir, "allAnnualPredictions.rds")
write_rds(data_pred, fname)
