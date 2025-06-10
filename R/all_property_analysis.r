data_test <- left_join(data_pred, bayes_estimates) |>
	rename(med_density = y)

data_test |>
	select(pred_density, med_density) |>
	filter(!is.na(med_density)) |>
	ggplot() +
	aes(x = med_density, y = pred_density) +
	geom_point() +
	geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 1) +
	labs(
		x = "Observed",
		y = "Predicted"
	) +
	theme_bw()

properties <- data_test |>
	pull(propertyID) |>
	unique()

all_slopes <- tibble()
for (j in seq_along(properties)) {
	tmp <- data_test |>
		filter(propertyID == properties[j]) |>
		arrange(year)

	if (nrow(tmp) == 1) {
		next
	}

	area <- tmp$property_area_km2[1]

	m <- lm(pred_density ~ year, data = tmp)
	gpr <- round(summary(m)$coefficients[2], 3)

	sp <- tibble(
		slope = gpr,
		st_name = tmp$st_name[1],
		propertyID = tmp$propertyID[1],
		N1 = round(tmp$pred_density[1] * area),
		county_code = tmp$county_code[1]
	)

	all_slopes <- bind_rows(all_slopes, sp)
}

good_states <- all_slopes |>
	count(st_name) |>
	filter(n >= 10) |>
	pull(st_name)

meds <- all_slopes |>
	filter(st_name %in% good_states) |>
	group_by(st_name) |>
	summarise(med = quantile(slope, 0.5)) |>
	ungroup()

all_slopes |>
	filter(st_name %in% good_states) |>
	group_by(st_name) |>
	summarise(
		low = quantile(slope, 0.05),
		q1 = quantile(slope, 0.25),
		med = quantile(slope, 0.5),
		q3 = quantile(slope, 0.75),
		high = quantile(slope, 0.95)
	) |>
	ggplot() +
	aes(x = med, y = reorder(st_name, med)) +
	geom_linerange(aes(xmin = low, xmax = high), linewidth = 0.5) +
	geom_linerange(aes(xmin = q1, xmax = q3), linewidth = 2.5) +
	geom_point(size = 5) +
	geom_point(size = 3, color = "white") +
	geom_vline(xintercept = 0, linetype = "dashed") +
	labs(
		x = element_blank(),
		y = "Yearly change in wild pig density",
		fill = "Management"
	) +
	theme_bw()

ggsave("plots/ML/annualChangeStateBoxWhisker.jpeg")


all_slopes |>
	filter(st_name %in% good_states) |>
	left_join(meds) |>
	ggplot() +
	aes(x = reorder(st_name, med), y = slope) +
	geom_boxplot(width = 0.5, alpha = 0.3, outliers = FALSE) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	labs(
		x = element_blank(),
		y = "Yearly change in wild pig density",
		fill = "Management"
	) +
	theme_bw() +
	theme(
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14)
	) +
	coord_flip()


property_cols <- rep(2, length(properties))
names(property_cols) <- properties

plot_sample_size <- function(df) {
	df |>
		select(propertyID, st_name, year) |>
		group_by(st_name, year) |>
		count() |>
		bind_rows(legend) |>
		ggplot() +
		aes(x = year, y = n, color = st_name) +
		geom_point() +
		geom_line() +
		labs(
			color = "State",
			x = "Year",
			y = "Number of properties"
		) +
		theme_bw() +
		theme(
			axis.text = element_text(size = 12),
			axis.title = element_text(size = 14)
		)
}

legend <- data_test |>
	select(st_name) |>
	distinct() |>
	mutate(n = NA)

s1 <- data_test |>
	filter(st_name == "TEXAS") |>
	plot_sample_size() +
	labs(title = "Texas", x = "")

s2 <- data_test |>
	filter(st_name != "TEXAS") |>
	plot_sample_size() +
	labs(title = "All other states")

ggarrange(
	s1,
	s2,
	nrow = 2,
	align = "hv",
	heights = c(1.5, 2),
	common.legend = TRUE,
	legend = "right"
)

ggsave("plots/ML/nPropertiesPerStateTimeSeries.jpeg")


x_breaks <- seq(2015, 2024, by = 3)

data_test |>
	ggplot() +
	aes(x = pred_density) +
	geom_density() +
	facet_wrap(~st_name) +
	theme_bw()

data_test |>
	select(propertyID, st_name, year, pred_density, all_take) |>
	group_by(st_name, year) |>
	summarise(`Median density` = median(pred_density), take = sum(all_take)) |>
	ungroup() |>
	mutate(`Log(total pigs removed + 1)` = log(take + 1)) |>
	pivot_longer(
		cols = c(`Log(total pigs removed + 1)`, `Median density`),
		names_to = "metric",
		values_to = "y"
	) |>
	ggplot() +
	aes(x = year, y = y, color = metric) +
	geom_point() +
	geom_line() +
	facet_wrap(~st_name) +
	labs(
		x = "Year",
		y = "Value",
		color = "Data"
	) +
	scale_x_continuous(breaks = x_breaks) +
	theme_bw() +
	theme(
		axis.text = element_text(size = 12),
		axis.text.x = element_text(angle = 45, vjust = 0.5),
		axis.title = element_text(size = 14)
	)

ggsave("plots/ML/densityTakeTimeSeriesByState.jpeg")

library(usmap)

df1 <- data_test |>
	select(county_code) |>
	unique() |>
	mutate(fips = as.character(county_code)) |>
	mutate(values = 0)

plot_usmap(
	data = df1,
	regions = "counties",
	exclude = c("AK", "HI")
) +
	scale_fill_continuous(
		low = "white",
		high = "red",
		name = "Counties with at least 1 density estimate",
		label = scales::comma
	) +
	theme(legend.position = "right")
