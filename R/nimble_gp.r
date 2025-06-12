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
		sigma2 <- sigma * sigma
		for (i in 1:n) {
			for (j in 1:n) {
				if (i == j) {
					result[i, j] <- sigma2 + 1e-6 # add a small value to the diagonal for numerical stability
				} else {
					result[i, j] <- pow(tau, 2) * exp(-rho * pow(dists[i, j], 2))
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
	rho_s ~ dunif(0, 5) # decay parameter for the GP autocorrelation function
	sigma_s ~ dunif(0, 10) # nugget
	tau_s ~ dunif(0, 10)

	# priors - time
	mu_t ~ dnorm(0, 0.01)
	rho_t ~ dunif(0, 5)
	sigma_t ~ dunif(0, 10)
	tau_t ~ dunif(0, 10)

	# priors - observation error
	tau_obs ~ dgamma(0.001, 0.001)

	# spatial GP
	alpha_s[1:N] <- mu_s * ones[1:N]
	cov_s[1:N, 1:N] <- expcov(dist_s[1:N, 1:N], rho_s, sigma_s, tau_s)
	s[1:N] ~ dmnorm(alpha_s[1:N], cov = cov_s[1:N, 1:N])

	# temporal GP
	alpha_t[1:N] <- mu_t * ones[1:N]
	cov_t[1:N, 1:N] <- expcov(dist_t[1:N, 1:N], rho_t, sigma_t, tau_t)
	t[1:N] ~ dmnorm(alpha_t[1:N], cov = cov_t[1:N, 1:N])

	# likelihood
	for (i in 1:N) {
		# combine spatial and temporal effects
		mu_x[i] <- max(0, s[i] + t[i]) # ensure positive values

		# observation model
		y[i] ~ dnorm(mu_x[i], tau = tau_obs)
	}
})
