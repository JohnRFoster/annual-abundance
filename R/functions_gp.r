library(nimble)

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
	# priors - space
	mu_s ~ dnorm(0, 0.01) # mean
	rho_s ~ dunif(0, 100) # decay parameter for the GP autocorrelation function
	sigma_s ~ dunif(0, 100) # nugget
	tau_cov_s ~ dunif(0, 100) # GP standard deviation parameter

	# priors - time
	mu_t ~ dnorm(0, 0.01)
	rho_t ~ dunif(0, 100)
	sigma_t ~ dunif(0, 100)
	tau_cov_t ~ dunif(0, 100)

	# priors - observation error
	tau_obs ~ dgamma(0.001, 0.001)

	# spatial GP
	alpha_s[1:N] <- mu_s * ones[1:N]
	cov_s[1:N, 1:N] <- expcov(dist_s[1:N, 1:N], rho_s, sigma_s, tau_cov_s)
	s[1:N] ~ dmnorm(alpha_s[1:N], cov = cov_s[1:N, 1:N])

	# temporal GP
	alpha_t[1:N] <- mu_t * ones[1:N]
	cov_t[1:N, 1:N] <- expcov(dist_t[1:N, 1:N], rho_t, sigma_t, tau_cov_t)
	t[1:N] ~ dmnorm(alpha_t[1:N], cov = cov_t[1:N, 1:N])

	# likelihood
	for (i in 1:N) {
		log(mu_x[i]) <- s[i] + t[i] # combine spatial and temporal effects
		y[i] ~ dnorm(mu_x[i], tau = tau_obs)
	}
})

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

create_spatial_distance_matrix <- function(df) {
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

	dist_matrix <- geosphere::distm(latlon)

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

	# normalize the matrix
	dist_matrix / max(dist_matrix)
}

create_year_distance_matrix <- function(df) {
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

	# Normalize the year matrix
	year_matrix / max(year_matrix)
}
