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
#'   `length` (structural length in cm, i.e., \eqn{V^{1/3}}).
#'   Optional additional column `weight` (wet mass in g).
#' @param reproduction A data frame with columns: `id`, `t_start`, `t_end`,
#'   `count` (number of offspring in the interval).
#'   For cumulative counts use [repro_to_intervals()] first.
#' @param survival A data frame with columns: `id`, `time`, `event`
#'   (1 = death observed, 0 = right-censored).
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
                      survival = NULL,
                      concentration = NULL,
                      f_food = 1.0) {

	if (is.null(growth) && is.null(reproduction) && is.null(survival)) {
		cli::cli_abort("At least one of {.arg growth}, {.arg reproduction}, or {.arg survival} must be provided.")
	}

	out <- list(
		growth       = NULL,
		reproduction = NULL,
		survival     = NULL,
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

	# --- Survival data ---
	if (!is.null(survival)) {
		survival <- validate_survival(survival)
		out$survival <- survival
		out$endpoints <- c(out$endpoints, "survival")
	}

	# Determine unique individuals across all endpoints
	all_ids <- character(0)
	if (!is.null(growth))       all_ids <- union(all_ids, as.character(growth$id))
	if (!is.null(reproduction)) all_ids <- union(all_ids, as.character(reproduction$id))
	if (!is.null(survival))     all_ids <- union(all_ids, as.character(survival$id))

	out$ids   <- sort(all_ids)
	out$n_ind <- length(out$ids)

	structure(out, class = "bdeb_data")
}

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
	if (!is.null(x$survival)) {
		n_events <- sum(x$survival$event == 1)
		cli::cli_alert("Survival: {n_events} deaths observed")
	}
	if (!is.null(x$concentration)) {
		cli::cli_alert("Concentration groups: {length(unique(x$concentration))}")
	}
	invisible(x)
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

validate_survival <- function(df) {
	required <- c("id", "time", "event")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		cli::cli_abort("Survival data missing columns: {.field {missing}}")
	}
	if (!all(df$event %in% c(0, 1))) {
		cli::cli_abort("Survival {.field event} must be 0 (censored) or 1 (death).")
	}
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
repro_to_intervals <- function(df) {
	required <- c("id", "time", "cumulative")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		cli::cli_abort("Data missing columns: {.field {missing}}")
	}
	df <- df[order(df$id, df$time), ]

	out_list <- lapply(split(df, df$id), function(d) {
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

#' Build Stan Data List for Individual Growth
#' @param data A `bdeb_data` object.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_individual <- function(data, priors) {
	g <- data$growth
	ids <- unique(g$id)

	if (length(ids) > 1) {
		cli::cli_abort("Individual model expects single individual. Use hierarchical for multiple.")
	}

	g <- g[order(g$time), ]

	stan_data <- list(
		N_obs  = nrow(g),
		t_obs  = g$time,
		L_obs  = g$length,
		f_food = data$f_food
	)

	# Append prior hyperparameters
	c(stan_data, prior_to_stan_data(priors))
}

#' Build Stan Data List for Hierarchical Growth
#' @param data A `bdeb_data` object.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_hierarchical <- function(data, priors) {
	g <- data$growth
	ids <- unique(g$id)
	n_ind <- length(ids)

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

	c(stan_data, prior_to_stan_data_hierarchical(priors))
}

#' Build Stan Data List for Growth + Reproduction
#' @param data A `bdeb_data` object with growth and reproduction.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_growth_repro <- function(data, priors) {
	g <- data$growth
	r <- data$reproduction

	# All unique times needed for ODE solving
	all_times <- sort(unique(c(g$time, r$t_start, r$t_end)))
	all_times <- all_times[all_times > 0]

	# Index mappings
	idx_L <- match(g$time, all_times)
	idx_R_start <- match(r$t_start, all_times)
	idx_R_end   <- match(r$t_end, all_times)

	# Handle t=0 observations by shifting to first positive time
	if (any(g$time == 0)) {
		t0_obs <- which(g$time == 0)
		all_times <- sort(unique(c(1e-3, all_times)))
		idx_L <- match(g$time, all_times)
		idx_L[t0_obs] <- 1
	}

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

	c(stan_data, prior_to_stan_data_growth_repro(priors))
}

#' Build Stan Data List for DEBtox
#' @param data A `bdeb_data` object with concentration info.
#' @param priors A list of `bdeb_prior` objects.
#' @return Named list suitable for Stan.
#' @keywords internal
build_stan_data_debtox <- function(data, priors) {
	g <- data$growth
	conc <- data$concentration

	if (is.null(conc)) {
		cli::cli_abort("DEBtox model requires {.arg concentration} in bdeb_data().")
	}

	# Group by concentration
	if (is.data.frame(conc)) {
		g <- merge(g, conc, by = "id")
	} else {
		g$concentration <- conc[match(g$id, names(conc))]
	}

	groups <- sort(unique(g$concentration))
	n_groups <- length(groups)

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

	# Reproduction (optional)
	has_repro <- as.integer(!is.null(data$reproduction))
	max_N_R <- 0L
	N_R_vec <- rep(0L, n_groups)
	R_counts_mat <- matrix(0L, nrow = n_groups, ncol = 1)
	idx_R_end_mat <- matrix(1L, nrow = n_groups, ncol = 1)

	if (has_repro == 1 && !is.null(data$reproduction)) {
		r <- data$reproduction
		# Map reproduction to groups (simplified)
		r_by_group <- split(r, r$id)  # placeholder — needs concentration mapping
		# For now, set empty
		has_repro <- 0L
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
		idx_R_end   = idx_R_end_mat
	)

	c(stan_data, prior_to_stan_data_debtox(priors))
}
