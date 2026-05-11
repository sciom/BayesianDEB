suppressPackageStartupMessages({
	library(BayesianDEB)
	library(posterior)
})

set.seed(42)

# Replicate 7.1 fit (Neuhauser 1980 growth) exactly as 03_case_studies.R
curves_path <- file.path("data", "curves.txt")
curves <- read.delim(curves_path)
neuh <- curves[curves$CURVE_ID == "gr0226", c("time", "bm")]
names(neuh) <- c("time", "wet_mass_mg")
d_V <- 1.05
neuh$length <- ((neuh$wet_mass_mg / 1000) / d_V)^(1 / 3)
eisenia_real <- data.frame(id = 1, time = neuh$time, length = neuh$length)

dat_real <- bdeb_data(growth = eisenia_real, f_food = 1.0)
mod_real <- bdeb_model(dat_real, type = "individual",
	priors = list(
		p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
		p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
		kappa   = prior_beta(a = 3, b = 2),
		v       = prior_lognormal(mu = -1.5, sigma = 0.5),
		E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
		sigma_L = prior_halfnormal(sigma = 0.05)),
	temperature = list(T_obs = 298.15, T_ref = 293.15, T_A = 8000))

fit_real <- bdeb_fit(mod_real,
	chains = 2, iter_warmup = 300, iter_sampling = 300,
	adapt_delta = 0.95, max_treedepth = 12, seed = 42, refresh = 0)

draws <- fit_real$fit$draws()
df <- posterior::as_draws_df(draws)
cat("\n--- nrow:", nrow(df), "chains:",
	paste(sort(unique(df$.chain)), collapse = ","), "---\n")

# NaN count per variable per chain
cat("\n--- NaN count by variable × chain ---\n")
vars <- setdiff(names(df), c(".chain", ".iteration", ".draw"))
res <- list()
for (v in vars) {
	for (ch in sort(unique(df$.chain))) {
		col <- df[[v]][df$.chain == ch]
		n_na <- sum(is.na(col))
		if (n_na > 0) {
			res[[length(res) + 1L]] <- data.frame(var = v, chain = ch, n_na = n_na)
		}
	}
}
if (length(res)) {
	r <- do.call(rbind, res)
	print(r, row.names = FALSE)
} else {
	cat("No NaN values in any variable.\n")
}

cat("\n--- Core pars (those passed to plot trace) ---\n")
core <- c("p_Am", "p_M", "kappa", "sigma_L")
for (v in intersect(core, vars)) {
	for (ch in sort(unique(df$.chain))) {
		x <- df[[v]][df$.chain == ch]
		cat(sprintf("%-10s chain %d: NA=%d  finite=%d  range=[%g, %g]\n",
			v, ch, sum(is.na(x)), sum(is.finite(x)),
			suppressWarnings(min(x, na.rm = TRUE)),
			suppressWarnings(max(x, na.rm = TRUE))))
	}
}

cat("\n--- Try plot(fit_real, type='trace', pars=c('p_Am','p_M','kappa','sigma_L')) ---\n")
out <- tryCatch({
	p <- plot(fit_real, type = "trace",
	          pars = c("p_Am", "p_M", "kappa", "sigma_L"))
	cat("plot() succeeded.\n")
	"OK"
}, error = function(e) {
	cat("plot() ERROR:\n")
	print(e)
	"FAIL"
})
