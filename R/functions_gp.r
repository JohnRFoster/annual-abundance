make_all_prop_years_gp <- function(df) {
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
			select(propertyID, st_name, Lat, Long)

		tmpp <- tibble(
			propertyID = prop_vec[i],
			year = min_year:max_year
		) |>
			left_join(prop_info, by = "propertyID")

		all_timesteps <- bind_rows(all_timesteps, tmpp)
	}

	all_timesteps
}

create_spatial_distance_matrix <- function(df, normalize) {
	if (!all(c("Long", "Lat") %in% names(df))) {
		stop("Data frame must contain 'Long' and 'Lat' columns.")
	}

	# Ensure the coordinates are numeric
	df$Long <- as.numeric(df$Long)
	df$Lat <- as.numeric(df$Lat)

	if (any(is.na(df$Long)) || any(is.na(df$Lat))) {
		stop("Coordinates contain NA values.")
	}

	# Calculate the distance matrix using geosphere
	latlon <- df |>
		select(Long, Lat) |>
		as.matrix()

	dist_matrix <- geosphere::distm(latlon) / 1000 # Convert to kilometers

	assertthat::are_equal(
		nrow(dist_matrix),
		nrow(df),
		msg = "Distance matrix dimensions do not match the number of rows in the data frame."
	)
	assertthat::are_equal(
		ncol(dist_matrix),
		nrow(df),
		msg = "Distance matrix dimensions do not match the number of rows in the data frame."
	)

	if (normalize) {
		dist_matrix / max(dist_matrix)
	} else {
		dist_matrix
	}
}

create_year_distance_matrix <- function(df, normalize) {
	if (!"year" %in% names(df)) {
		stop("Data frame must contain a 'year' column.")
	}

	# Ensure the year column is numeric
	df$year <- as.numeric(df$year)

	if (any(is.na(df$year))) {
		stop("Year column contains NA values.")
	}

	# Create a matrix of years
	x <- matrix(df$year, ncol = 1)
	year_matrix <- as.matrix(dist(x, method = "maximum"))

	assertthat::are_equal(
		nrow(year_matrix),
		nrow(df),
		msg = "Temporal distance matrix dimensions do not match the number of rows in the data frame."
	)
	assertthat::are_equal(
		ncol(year_matrix),
		nrow(df),
		msg = "Temporal distance matrix dimensions do not match the number of rows in the data frame."
	)

	if (normalize) {
		year_matrix / max(year_matrix)
	} else {
		year_matrix
	}
}
