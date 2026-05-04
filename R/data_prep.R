#' Prepare Data for BDEB Models
#'
#' Converts long-format data frames into the structured list required by
#' the BayesianDEB Stan programs.  Growth observations are matched to the
#' DEB state variable \eqn{L = V^{1/3}} (structural length); reproduction
#' records are interval counts of offspring over \eqn{[t_{\mathrm{start}},
#' t_{\mathrm{end}})}.  The function validates column names, rejects
#' negative times/lengths, sorts by individual and time, and (for
#' hierarchical models) pads ragged observation vectors into matrices with
#' `NaN` fill, as required by Stan's fixed-size array declarations.
#'
#' @param growth A data frame with columns: `id`, `time` (days),
#'   `length` (**structural** length in cm, i.e., \eqn{L = V^{1/3}}).
#'   If your data are **physical (observed) lengths**, multiply by the
#'   shape coefficient: \eqn{L = \delta_M \times L_w}.  See the
#'   package vignette for guidance on \eqn{\delta_M}.
#' @param reproduction A data frame with columns: `id`, `t_start`, `t_end`,
#'   `count` (number of offspring in the interval).
#'   For cumulative counts use [repro_to_intervals()] first.
#' @param concentration Optional named numeric vector or data frame mapping
#'   individual/group id to external toxicant concentration \eqn{C_w}
#'   (for DEBtox models; see [bdeb_tox()]).
#' @param f_food Scaled functional response \eqn{f \in [0,1]}, the ratio
#'   of actual to maximum ingestion rate (Kooijman, 2010, Eq. 2.3).
#'   Default 1 (ad libitum feeding).
#' @return A `bdeb_data` object (S3 list) ready for [bdeb_model()].
#'
#' @references
#' Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
#' Organisation*. 3rd edition. Cambridge University Press.
#' \doi{10.1017/CBO9780511805400}
#' @export
#' @examples
#' # Simple growth data
#' df <- data.frame(
#'   id = rep(1, 10),
#'   time = seq(0, 45, by = 5),
#'   length = c(0.1, 0.15, 0.22, 0.30, 0.38, 0.44, 0.49, 0.52, 0.54, 0.55)
#' )
#' dat <- bdeb_data(growth = df)
bdeb_data <- function(growth = NULL,
                      reproduction = NULL,
                      concentration = NULL,
                      f_food = 1.0) {

	if (is.null(growth) && is.null(reproduction)) {
		cli::cli_abort("At least one of {.arg growth} or {.arg reproduction} must be provided.")
	}
	if (!is.null(growth) && !is.data.frame(growth)) {
		cli::cli_abort("{.arg growth} must be a data frame or NULL.")
	}
	if (!is.null(reproduction) && !is.data.frame(reproduction)) {
		cli::cli_abort("{.arg reproduction} must be a data frame or NULL.")
	}
	if (!is.numeric(f_food) || length(f_food) != 1 ||
	    !is.finite(f_food) || f_food < 0 || f_food > 1) {
		cli::cli_abort("{.arg f_food} must be a finite scalar in [0, 1].")
	}

	out <- list(
		growth       = NULL,
		reproduction = NULL,
		concentration = concentration,
		f_food       = f_food,
		n_ind        = 0L,
		ids          = character(0),
		endpoints    = character(0)
	)

	# --- Growth data ---
	if (!is.null(growth)) {
		growth <- validate_growth(growth)
		out$growth <- growth
		out$endpoints <- c(out$endpoints, "growth")
	}

	# --- Reproduction data ---
	if (!is.null(reproduction)) {
		reproduction <- validate_repro(reproduction)
		out$reproduction <- reproduction
		out$endpoints <- c(out$endpoints, "reproduction")
	}

	# Determine unique individuals across all endpoints
	all_ids <- character(0)
	if (!is.null(growth))       all_ids <- union(all_ids, as.character(growth$id))
	if (!is.null(reproduction)) all_ids <- union(all_ids, as.character(reproduction$id))

	out$ids   <- sort(all_ids)
	out$n_ind <- length(out$ids)

	structure(out, class = "bdeb_data")
}

#' Print a BDEB Data Object
#'
#' @param x A [bdeb_data()] object.
#' @param ... Ignored.
#' @return The input object, invisibly.
#' @export
print.bdeb_data <- function(x, ...) {
	cli::cli_h2("BDEB Data")
	cli::cli_alert_info("Individuals: {x$n_ind}")
	cli::cli_alert_info("Endpoints: {paste(x$endpoints, collapse = ', ')}")
	cli::cli_alert_info("Functional response (f): {x$f_food}")
	if (!is.null(x$growth)) {
		n_obs <- nrow(x$growth)
		t_range <- range(x$growth$time)
		cli::cli_alert("Growth: {n_obs} observations, t = [{t_range[1]}, {t_range[2]}]")
	}
	if (!is.null(x$reproduction)) {
		n_obs <- nrow(x$reproduction)
		cli::cli_alert("Reproduction: {n_obs} interval records")
	}
	if (!is.null(x$concentration)) {
		cli::cli_alert("Concentration groups: {length(unique(x$concentration))}")
	}
	invisible(x)
}

#' Summary of a BDEB Data Object
#'
#' Returns a compact list of summary statistics describing a
#' [bdeb_data()] object: number of observations, number of unique
#' individuals, time range, presence of growth/reproduction endpoints,
#' the functional response, and (for DEBtox data) the unique
#' concentration levels.
#'
#' @param object A [bdeb_data()] object.
#' @param ... Ignored.
#' @return An object of class `summary.bdeb_data` (a list).
#' @export
#' @examples
#' df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
#' dat <- bdeb_data(growth = df)
#' summary(dat)
summary.bdeb_data <- function(object, ...) {
	if (!inherits(object, "bdeb_data")) {
		cli::cli_abort("{.arg object} must be a {.cls bdeb_data} object.")
	}

	n_growth <- if (!is.null(object$growth)) nrow(object$growth) else 0L
	n_repro  <- if (!is.null(object$reproduction)) nrow(object$reproduction) else 0L

	t_range <- NULL
	if (!is.null(object$growth)) {
		t_range <- range(object$growth$time, na.rm = TRUE)
	} else if (!is.null(object$reproduction)) {
		t_range <- range(c(object$reproduction$t_start,
		                   object$reproduction$t_end), na.rm = TRUE)
	}

	out <- list(
		n_obs_growth     = n_growth,
		n_obs_repro      = n_repro,
		n_individuals    = object$n_ind,
		ids              = object$ids,
		time_range       = t_range,
		has_growth       = !is.null(object$growth),
		has_reproduction = !is.null(object$reproduction),
		f_food           = object$f_food,
		endpoints        = object$endpoints
	)

	if (!is.null(object$concentration)) {
		conc <- if (is.data.frame(object$concentration)) {
			object$concentration$concentration
		} else {
			as.numeric(object$concentration)
		}
		out$conc_levels <- sort(unique(conc))
	}

	structure(out, class = "summary.bdeb_data")
}

#' @rdname summary.bdeb_data
#' @param x A `summary.bdeb_data` object.
#' @export
print.summary.bdeb_data <- function(x, ...) {
	cli::cli_h2("BDEB Data Summary")
	cli::cli_li("Individuals: {.val {x$n_individuals}}")
	if (x$has_growth) {
		cli::cli_li("Growth observations: {.val {x$n_obs_growth}}")
	}
	if (x$has_reproduction) {
		cli::cli_li("Reproduction records: {.val {x$n_obs_repro}}")
	}
	if (!is.null(x$time_range)) {
		cli::cli_li(
			"Time range: [{format(x$time_range[1], digits = 3)}, {format(x$time_range[2], digits = 3)}]"
		)
	}
	cli::cli_li("Functional response (f): {.val {x$f_food}}")
	if (!is.null(x$conc_levels)) {
		cli::cli_li(
			"Concentration levels: {paste(x$conc_levels, collapse = ', ')}"
		)
	}
	invisible(x)
}

#' Plot a BDEB Data Object
#'
#' Visualises the contents of a [bdeb_data()] object.  Growth data are
#' shown as observed length versus time, with one trace per individual.
#' Reproduction data are shown as interval counts versus interval
#' midpoint.  When concentration information is available (DEBtox
#' setups) individuals are coloured by group.
#'
#' @param x A [bdeb_data()] object.
#' @param endpoint Which endpoint to plot: `"growth"` or
#'   `"reproduction"`.  Default `NULL`, which prefers growth when both
#'   are present.
#' @param ... Ignored.
#' @return A [ggplot2::ggplot] object.
#' @export
#' @examples
#' df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
#' dat <- bdeb_data(growth = df)
#' plot(dat)
plot.bdeb_data <- function(x, endpoint = NULL, ...) {
	if (!inherits(x, "bdeb_data")) {
		cli::cli_abort("{.arg x} must be a {.cls bdeb_data} object.")
	}

	if (is.null(endpoint)) {
		endpoint <- if (!is.null(x$growth)) "growth" else "reproduction"
	}
	endpoint <- match.arg(endpoint, c("growth", "reproduction"))

	if (endpoint == "growth") {
		if (is.null(x$growth)) {
			cli::cli_abort("No growth data to plot.")
		}
		df <- x$growth
		df$id <- factor(df$id)

		if (!is.null(x$concentration) && !("concentration" %in% names(df))) {
			conc <- x$concentration
			if (is.data.frame(conc)) {
				df <- merge(df, conc, by = "id", all.x = TRUE)
			} else {
				df$concentration <- conc[match(as.character(df$id),
				                               names(conc))]
			}
		}

		p <- ggplot2::ggplot(df,
		                    ggplot2::aes(x = .data$time,
		                                 y = .data$length,
		                                 group = .data$id))

		if ("concentration" %in% names(df) &&
		    any(!is.na(df$concentration))) {
			p <- p + ggplot2::aes(colour = factor(.data$concentration)) +
				ggplot2::labs(colour = "C_w")
		}

		p +
			ggplot2::geom_point() +
			ggplot2::geom_line(alpha = 0.4) +
			ggplot2::theme_bw() +
			ggplot2::labs(
				title = "BDEB growth data",
				x = "Time", y = "Structural length"
			)
	} else {
		if (is.null(x$reproduction)) {
			cli::cli_abort("No reproduction data to plot.")
		}
		df <- x$reproduction
		df$id <- factor(df$id)
		df$t_mid <- 0.5 * (df$t_start + df$t_end)

		ggplot2::ggplot(df,
		                ggplot2::aes(x = .data$t_mid, y = .data$count,
		                             group = .data$id)) +
			ggplot2::geom_point() +
			ggplot2::geom_line(alpha = 0.4) +
			ggplot2::theme_bw() +
			ggplot2::labs(
				title = "BDEB reproduction data",
				x = "Interval midpoint", y = "Offspring count"
			)
	}
}

# --- Internal validation helpers ---

validate_growth <- function(df) {
	required <- c("id", "time", "length")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		cli::cli_abort("Growth data missing columns: {.field {missing}}")
	}
	df <- df[order(df$id, df$time), ]
	if (any(df$time < 0))  cli::cli_abort("Growth times must be non-negative.")
	if (any(df$length < 0, na.rm = TRUE)) cli::cli_abort("Growth lengths must be non-negative.")

	# Warn if lengths look like physical (observed) rather than structural
	max_L <- max(df$length, na.rm = TRUE)
	if (max_L > 10) {
		cli::cli_warn(c(
			"!" = "Maximum length is {round(max_L, 1)} cm, which is unusually large for structural length.",
			"i" = "BayesianDEB expects {.emph structural} length L = V^(1/3), not physical length.",
			"i" = "If these are physical lengths, convert: L_structural = delta_M * L_physical.",
			"i" = "Typical delta_M values: 0.1-0.5 (see AmP entry for your species).",
			"i" = "Suppress this warning if your data are truly structural lengths > 10 cm."
		))
	}
	df
}

validate_repro <- function(df) {
	required <- c("id", "t_start", "t_end", "count")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		cli::cli_abort("Reproduction data missing columns: {.field {missing}}")
	}
	df <- df[order(df$id, df$t_start), ]
	if (any(df$t_end <= df$t_start)) {
		cli::cli_abort("Reproduction {.field t_end} must be > {.field t_start}.")
	}
	if (any(df$count < 0)) cli::cli_abort("Reproduction counts must be non-negative.")
	df
}

#' Convert Cumulative Reproduction to Intervals
#'
#' Many ecotoxicological protocols (e.g., ISO 11267 for *Folsomia candida*,
#' OECD 222 for *Eisenia fetida*) report cumulative offspring counts at
#' successive observation times.  The DEB reproduction buffer model, however,
#' requires interval counts \eqn{\Delta R = R(t_{\mathrm{end}}) -
#' R(t_{\mathrm{start}})} so that the negative-binomial likelihood can be
#' applied to each counting period.  This function computes the
#' first-difference per individual.
#'
#' @param df Data frame with columns: `id`, `time`, `cumulative`.
#' @return Data frame with columns: `id`, `t_start`, `t_end`, `count`.
#' @export
#' @examples
#' cumul <- data.frame(
#'   id = rep(1, 5),
#'   time = c(0, 7, 14, 21, 28),
#'   cumulative = c(0, 10, 30, 60, 100)
#' )
#' repro_to_intervals(cumul)
repro_to_intervals <- function(df) {
	required <- c("id", "time", "cumulative")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		cli::cli_abort("Data missing columns: {.field {missing}}")
	}
	df <- df[order(df$id, df$time), ]

	by_id <- split(df, df$id)
	single_obs <- names(which(vapply(by_id, nrow, integer(1)) < 2))
	if (length(single_obs) > 0) {
		cli::cli_warn(
			"Dropping {length(single_obs)} individual(s) with < 2 observations: {.val {single_obs}}. Need at least 2 time points to compute intervals."
		)
	}

	out_list <- lapply(by_id, function(d) {
		n <- nrow(d)
		if (n < 2) return(NULL)
		data.frame(
			id      = d$id[-1],
			t_start = d$time[-n],
			t_end   = d$time[-1],
			count   = diff(d$cumulative)
		)
	})
	do.call(rbind, out_list)
}

#' Encode observation model flags for Stan
#' @param observation Named list of bdeb_obs objects.
#' @return Named list with obs_growth, obs_nu, obs_repro.
#' @keywords internal
observation_to_stan_data <- function(observation) {
	if (is.null(observation)) {
		observation <- list(growth = obs_normal(), reproduction = obs_negbinom())
	}
	# Growth: normal=1, lognormal=2, student_t=3
	g_family <- observation$growth$family
	obs_growth <- switch(g_family,
		normal    = 1L,
		lognormal = 2L,
		student_t = 3L,
		1L  # fallback
	)
	obs_nu <- if (g_family == "student_t") observation$growth$nu else 5.0

	# Reproduction: negbinom=1, poisson=2
	r_family <- observation$reproduction$family
	obs_repro <- switch(r_family,
		negbinom = 1L,
		poisson  = 2L,
		1L  # fallback
	)

	list(obs_growth = obs_growth, obs_nu = obs_nu, obs_repro = obs_repro)
}

#' Encode temperature correction for Stan
#' @param temperature NULL or list with T_obs, T_ref, T_A.
#' @return Named list with has_temperature, T_obs, T_ref, T_A.
#' @keywords internal
temperature_to_stan_data <- function(temperature) {
	if (is.null(temperature)) {
		list(has_temperature = 0L, T_obs = 293.15, T_ref = 293.15, T_A = 0.0)
	} else {
		list(
			has_temperature = 1L,
			T_obs = temperature$T_obs,
			T_ref = temperature$T_ref,
			T_A   = temperature$T_A
		)
	}
}

#' Build Stan Data List for Individual Growth
#' @param data A `bdeb_data` object.
#' @param priors A list of `bdeb_prior` objects.
#' @param temperature NULL or list with T_obs, T_ref, T_A.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_individual <- function(data, priors, temperature = NULL,
                                      observation = NULL) {
	g <- data$growth
	ids <- unique(g$id)

	if (length(ids) > 1) {
		cli::cli_abort("Individual model expects single individual. Use hierarchical for multiple.")
	}

	g <- g[order(g$time), ]

	# ODE solver requires all t_obs > 0 (initial time is t0=0).
	# Replace t=0 with small epsilon; corresponding L_obs constrains L0.
	if (any(g$time == 0)) {
		g$time[g$time == 0] <- 1e-3
	}

	stan_data <- list(
		N_obs  = nrow(g),
		t_obs  = g$time,
		L_obs  = g$length,
		f_food = data$f_food
	)

	c(stan_data, observation_to_stan_data(observation),
	  temperature_to_stan_data(temperature), prior_to_stan_data(priors))
}

#' Build Stan Data List for Hierarchical Growth
#' @param data A `bdeb_data` object.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_hierarchical <- function(data, priors, temperature = NULL,
                                        observation = NULL) {
	g <- data$growth
	ids <- unique(g$id)
	n_ind <- length(ids)

	# ODE solver requires all t_obs > 0 (initial time is t0=0).
	if (any(g$time == 0)) {
		g$time[g$time == 0] <- 1e-3
	}

	# Split by individual
	by_id <- split(g, g$id)
	n_obs <- vapply(by_id, nrow, integer(1))
	max_n_obs <- max(n_obs)

	# Pad into matrices
	t_mat <- matrix(0, nrow = n_ind, ncol = max_n_obs)
	L_mat <- matrix(NaN, nrow = n_ind, ncol = max_n_obs)

	for (j in seq_along(by_id)) {
		d <- by_id[[j]]
		d <- d[order(d$time), ]
		ni <- nrow(d)
		t_mat[j, 1:ni] <- d$time
		L_mat[j, 1:ni] <- d$length
	}

	stan_data <- list(
		N_ind     = n_ind,
		max_N_obs = max_n_obs,
		N_obs     = as.array(n_obs),
		t_obs     = t_mat,
		L_obs     = L_mat,
		f_food    = data$f_food
	)

	c(stan_data, observation_to_stan_data(observation),
	  temperature_to_stan_data(temperature),
	  prior_to_stan_data_hierarchical(priors))
}

#' Build Stan Data List for Growth + Reproduction
#' @param data A `bdeb_data` object with growth and reproduction.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_growth_repro <- function(data, priors, temperature = NULL,
                                        observation = NULL) {
	g <- data$growth
	r <- data$reproduction

	# ODE solver requires all t_obs > 0 (initial time is t0=0).
	if (any(g$time == 0)) g$time[g$time == 0] <- 1e-3
	if (any(r$t_start == 0)) r$t_start[r$t_start == 0] <- 1e-3
	if (any(r$t_end == 0)) r$t_end[r$t_end == 0] <- 1e-3

	# All unique times needed for ODE solving
	all_times <- sort(unique(c(g$time, r$t_start, r$t_end)))

	# Index mappings
	idx_L <- match(g$time, all_times)
	idx_R_start <- match(r$t_start, all_times)
	idx_R_end   <- match(r$t_end, all_times)

	stan_data <- list(
		N_L         = nrow(g),
		t_L         = g$time,
		L_obs       = g$length,
		N_R         = nrow(r),
		t_R_start   = r$t_start,
		t_R_end     = r$t_end,
		R_counts    = as.integer(r$count),
		f_food      = data$f_food,
		N_times     = length(all_times),
		t_all       = all_times,
		idx_L       = as.array(idx_L),
		idx_R_start = as.array(idx_R_start),
		idx_R_end   = as.array(idx_R_end)
	)

	c(stan_data, observation_to_stan_data(observation),
	  temperature_to_stan_data(temperature),
	  prior_to_stan_data_growth_repro(priors))
}

#' Build Stan Data List for DEBtox
#' @param data A `bdeb_data` object with concentration info.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_debtox <- function(data, priors, temperature = NULL,
                                  observation = NULL) {
	g <- data$growth
	conc <- data$concentration

	if (is.null(conc)) {
		cli::cli_abort("DEBtox model requires {.arg concentration} in bdeb_data().")
	}

	# Assign concentration to each growth record
	if ("concentration" %in% names(g)) {
		# Dataset already has a concentration column (e.g. debtox_growth)
	} else if (is.data.frame(conc)) {
		g <- merge(g, conc, by = "id")
	} else {
		# Named vector: names are IDs, values are concentrations
		g$concentration <- conc[match(as.character(g$id), names(conc))]
	}

	groups <- sort(unique(g$concentration))
	n_groups <- length(groups)

	# Check for multiple individuals per concentration group
	# DEBtox fits ONE ODE per group — individual-level data must be
	# aggregated to group means per time point.
	by_group_raw <- split(g, g$concentration)
	has_multi <- vapply(by_group_raw, function(d) {
		length(unique(d$id)) > 1
	}, logical(1))

	if (any(has_multi)) {
		n_multi <- sum(has_multi)
		cli::cli_warn(c(
			"!" = "{n_multi} concentration group(s) contain multiple individuals.",
			"i" = "DEBtox fits one ODE per concentration group, not per individual.",
			"i" = "Aggregating to group means per time point. For individual-level TKTD, a hierarchical DEBtox extension is needed (not yet implemented)."
		))
		# Aggregate: mean length per (concentration, time)
		g <- stats::aggregate(length ~ concentration + time, data = g, FUN = mean)
	}

	# ODE solver requires all t_obs > 0 (initial time is t0=0).
	if (any(g$time == 0)) {
		g$time[g$time == 0] <- 1e-3
	}

	by_group <- split(g, g$concentration)
	n_obs <- vapply(by_group, nrow, integer(1))
	max_n_obs <- max(n_obs)

	t_mat <- matrix(0, nrow = n_groups, ncol = max_n_obs)
	L_mat <- matrix(NaN, nrow = n_groups, ncol = max_n_obs)

	for (i in seq_along(by_group)) {
		d <- by_group[[i]]
		d <- d[order(d$time), ]
		ni <- nrow(d)
		t_mat[i, 1:ni] <- d$time
		L_mat[i, 1:ni] <- d$length
	}

	# Reproduction (optional): map repro records to concentration groups
	has_repro <- 0L
	max_N_R <- 1L
	N_R_vec <- rep(0L, n_groups)
	R_counts_mat <- matrix(0L, nrow = n_groups, ncol = 1)
	idx_R_start_mat <- matrix(1L, nrow = n_groups, ncol = 1)
	idx_R_end_mat <- matrix(1L, nrow = n_groups, ncol = 1)

	if (!is.null(data$reproduction)) {
		r <- data$reproduction

		# Match t=0 shift applied to growth times
		if (any(r$t_start == 0)) r$t_start[r$t_start == 0] <- 1e-3
		if (any(r$t_end == 0))   r$t_end[r$t_end == 0]     <- 1e-3

		# Assign concentration to each repro record
		if ("concentration" %in% names(r)) {
			# Already has concentration column
		} else if (is.data.frame(conc)) {
			r <- merge(r, conc, by = "id")
		} else {
			r$concentration <- conc[match(as.character(r$id), names(conc))]
		}

		# Drop records with unmapped concentration (NA)
		r <- r[!is.na(r$concentration), ]

		if (nrow(r) > 0) {
			has_repro <- 1L

			# For each group, find which ODE time points correspond to t_end
			# The ODE is solved at t_obs[g, 1:N_obs[g]] for growth.
			# Repro t_end must map to one of these time points.
			r_by_group <- split(r, r$concentration)

			N_R_list <- integer(n_groups)
			for (gi in seq_len(n_groups)) {
				grp_name <- as.character(groups[gi])
				if (grp_name %in% names(r_by_group)) {
					N_R_list[gi] <- nrow(r_by_group[[grp_name]])
				}
			}

			max_N_R <- max(max(N_R_list), 1L)
			N_R_vec <- N_R_list
			R_counts_mat <- matrix(0L, nrow = n_groups, ncol = max_N_R)
			idx_R_start_mat <- matrix(1L, nrow = n_groups, ncol = max_N_R)
			idx_R_end_mat <- matrix(1L, nrow = n_groups, ncol = max_N_R)

			for (gi in seq_len(n_groups)) {
				grp_name <- as.character(groups[gi])
				if (!(grp_name %in% names(r_by_group))) next
				rg <- r_by_group[[grp_name]]
				rg <- rg[order(rg$t_end), ]
				ni <- nrow(rg)

				# Time vector for this group
				t_g <- t_mat[gi, seq_len(n_obs[gi])]

				for (ri in seq_len(ni)) {
					R_counts_mat[gi, ri] <- as.integer(rg$count[ri])
					# Strict match: repro times must exist in growth times
					idx_s <- match(rg$t_start[ri], t_g)
					idx_e <- match(rg$t_end[ri], t_g)
					if (is.na(idx_s) || is.na(idx_e)) {
						cli::cli_abort(c(
							"x" = "DEBtox reproduction times must match growth observation times.",
							"i" = "Group {.val {grp_name}}: t_start={rg$t_start[ri]} or t_end={rg$t_end[ri]} not found in growth times.",
							"i" = "Available times: {paste(t_g, collapse=', ')}"
						))
					}
					idx_R_start_mat[gi, ri] <- idx_s
					idx_R_end_mat[gi, ri] <- idx_e
				}
			}
		}
	}

	stan_data <- list(
		N_groups    = n_groups,
		C_w         = groups,
		max_N_obs   = max_n_obs,
		N_obs       = as.array(n_obs),
		t_obs       = t_mat,
		L_obs       = L_mat,
		f_food      = data$f_food,
		has_repro   = has_repro,
		N_R         = as.array(N_R_vec),
		max_N_R     = max(max_N_R, 1L),
		R_counts    = R_counts_mat,
		idx_R_start = idx_R_start_mat,
		idx_R_end   = idx_R_end_mat
	)

	c(stan_data, observation_to_stan_data(observation),
	  temperature_to_stan_data(temperature),
	  prior_to_stan_data_debtox(priors))
}
