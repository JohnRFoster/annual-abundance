samples <- read_rds("out/gp_samples.rds")
time_join <- read_rds("out/gp_data.rds")

unobserved <- which(is.na(time_join$pred_density))

idx <- 1
node <- paste0("mu_x[", idx, "]")
hist(samples[[node]], breaks = 50, main = paste0("Posterior of ", node))
abline(v = time_join$pred_density[idx], col = "red", lwd = 2)

mu_x <- samples[, grep("mu_x", colnames(samples))]
med_dens <- apply(mu_x, 2, median)

summary(apply(mu_x, 2, min))

plot(
	time_join$pred_density[-unobserved],
	med_dens[-unobserved],
	xlab = "Observed (xgBoost) density",
	ylab = "Predicted (GP) density",
	main = "Posterior Median of Predicted Densities"
)
abline(0, 1, col = "red", lwd = 2)

cor(time_join$pred_density[-unobserved], med_dens[-unobserved])


node <- "rho_t"
hist(samples[[node]], breaks = 50, main = paste0("Posterior of ", node))
