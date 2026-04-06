# ===========================================================
# Mock bdeb_fit objects for testing downstream functions
# without requiring cmdstanr / CmdStan
# ===========================================================

#' Create a minimal mock of CmdStanMCMC draws for testing
#' @param n_draws Number of posterior draws
#' @param type One of "individual", "hierarchical", "debtox"
#' @return A mock bdeb_fit object
mock_bdeb_fit <- function(n_draws = 100, type = "individual") {

  n_chains <- 2L
  n_per_chain <- n_draws %/% n_chains

  # Base DEB parameters (realistic draws)
  set.seed(42)
  pars <- data.frame(
    p_Am    = rlnorm(n_draws, 1.5, 0.1),
    p_M     = rlnorm(n_draws, -0.7, 0.1),
    kappa   = rbeta(n_draws, 8, 3),
    v       = rlnorm(n_draws, -1.5, 0.1),
    E_G     = rlnorm(n_draws, 6.0, 0.1),
    E0      = rlnorm(n_draws, 0.0, 0.1),
    L0      = rlnorm(n_draws, -2.0, 0.1),
    sigma_L = abs(rnorm(n_draws, 0.015, 0.003))
  )

  if (type == "hierarchical") {
    pars$mu_log_p_Am    <- rnorm(n_draws, 1.5, 0.05)
    pars$sigma_log_p_Am <- abs(rnorm(n_draws, 0.1, 0.02))
    for (j in 1:3) {
      pars[[paste0("p_Am_ind[", j, "]")]] <- rlnorm(n_draws, 1.5, 0.1)
      pars[[paste0("L0[", j, "]")]] <- rlnorm(n_draws, -2.0, 0.05)
    }
    pars$p_Am_new <- rlnorm(n_draws, 1.5, 0.15)
    pars$p_Am <- NULL
    pars$L0 <- NULL  # scalar L0 replaced by L0[j]
  }

  if (type == "debtox") {
    pars$k_d <- rlnorm(n_draws, -1.0, 0.1)
    pars$z_w <- rlnorm(n_draws, 2.7, 0.1)
    pars$b_w <- rlnorm(n_draws, -5.5, 0.2)
    pars$EC50 <- pars$z_w + 0.5 / pars$b_w
    pars$NEC  <- pars$z_w
  }

  # L_hat (predicted lengths) — 10 time points
  n_obs <- 10
  for (i in seq_len(n_obs)) {
    pars[[paste0("L_hat[", i, "]")]] <- 0.1 + 0.04 * i + rnorm(n_draws, 0, 0.005)
  }

  # L_rep (replicated lengths for PPC)
  for (i in seq_len(n_obs)) {
    pars[[paste0("L_rep[", i, "]")]] <- pars[[paste0("L_hat[", i, "]")]] +
      rnorm(n_draws, 0, 0.015)
  }

  # log_lik
  for (i in seq_len(n_obs)) {
    pars[[paste0("log_lik[", i, "]")]] <- rnorm(n_draws, -1, 0.5)
  }

  pars$lp__ <- rnorm(n_draws, -50, 5)

  # Chain info
  pars$.chain     <- rep(seq_len(n_chains), each = n_per_chain)
  pars$.iteration <- rep(seq_len(n_per_chain), n_chains)
  pars$.draw      <- seq_len(n_draws)

  draws_obj <- posterior::as_draws_df(pars)

  # Mock CmdStanMCMC-like object
  mock_cmdstan <- list(
    draws = function(format = NULL, ...) {
      if (!is.null(format) && format == "draws_matrix") {
        return(posterior::as_draws_matrix(draws_obj))
      }
      draws_obj
    },
    diagnostic_summary = function(quiet = FALSE) {
      list(
        num_divergent     = rep(0L, n_chains),
        num_max_treedepth = rep(0L, n_chains),
        ebfmi             = rep(1.1, n_chains)
      )
    }
  )

  # Stan data that matches
  stan_data <- list(
    N_obs  = n_obs,
    t_obs  = seq(0, 63, length.out = n_obs),
    L_obs  = 0.1 + 0.04 * seq_len(n_obs) + rnorm(n_obs, 0, 0.01),
    f_food = 1.0
  )

  if (type == "hierarchical") {
    n_hier_ind <- 3L
    stan_data <- list(
      N_ind     = n_hier_ind,
      max_N_obs = n_obs,
      N_obs     = rep(as.integer(n_obs), n_hier_ind),
      t_obs     = matrix(rep(seq(0, 63, length.out = n_obs), n_hier_ind),
                         nrow = n_hier_ind, byrow = TRUE),
      L_obs     = matrix(0.1 + 0.04 * rep(seq_len(n_obs), n_hier_ind) +
                         rnorm(n_obs * n_hier_ind, 0, 0.01),
                         nrow = n_hier_ind, byrow = TRUE),
      f_food    = 1.0
    )
  }

  if (type == "debtox") {
    stan_data <- list(
      N_groups  = 4L,
      C_w       = c(0, 20, 80, 200),
      max_N_obs = 7L,
      N_obs     = rep(7L, 4),
      t_obs     = matrix(rep(seq(0, 42, by = 7), 4), nrow = 4, byrow = TRUE),
      L_obs     = matrix(rnorm(28, 0.3, 0.05), nrow = 4),
      f_food    = 1.0
    )
  }

  # bdeb_data mock
  mock_data <- structure(list(
    growth   = data.frame(id = 1, time = stan_data$t_obs[1:min(n_obs, length(stan_data$t_obs))],
                          length = stan_data$L_obs[1:min(n_obs, length(stan_data$L_obs))]),
    n_ind    = 1L,
    ids      = "1",
    endpoints = "growth",
    f_food   = 1.0,
    concentration = if (type == "debtox") c("1"=0,"2"=20,"3"=80,"4"=200) else NULL
  ), class = "bdeb_data")

  # bdeb_model mock
  mock_model <- structure(list(
    data            = mock_data,
    type            = type,
    priors          = prior_default(type),
    observation     = list(growth = obs_normal(), reproduction = obs_negbinom()),
    temperature     = NULL,
    stan_model_name = switch(type,
      individual = "bdeb_individual_growth",
      hierarchical = "bdeb_hierarchical_growth",
      debtox = "bdeb_debtox"
    ),
    stan_data       = stan_data
  ), class = "bdeb_model")

  structure(list(
    fit           = mock_cmdstan,
    model         = mock_model,
    stan_model    = NULL,
    chains        = n_chains,
    iter_warmup   = 500L,
    iter_sampling = as.integer(n_per_chain),
    adapt_delta   = 0.8
  ), class = "bdeb_fit")
}


#' Create a mock bdeb_fit with BAD diagnostics for testing
#' @param n_divergent Number of divergent transitions per chain
#' @param n_max_treedepth Number of max treedepth hits per chain
#' @param low_ebfmi If TRUE, set E-BFMI below 0.3
#' @param bimodal If TRUE, create bimodal posterior (bad Rhat)
mock_bdeb_fit_bad <- function(n_divergent = 50,
                              n_max_treedepth = 20,
                              low_ebfmi = TRUE,
                              bimodal = FALSE) {
  n_draws <- 200
  n_chains <- 2L
  n_per_chain <- n_draws %/% n_chains

  set.seed(99)

  if (bimodal) {
    # Two chains stuck at different modes => bad Rhat
    pars <- data.frame(
      p_Am    = c(rlnorm(n_per_chain, 1.0, 0.05),
                  rlnorm(n_per_chain, 2.5, 0.05)),
      p_M     = c(rlnorm(n_per_chain, -0.5, 0.05),
                  rlnorm(n_per_chain, 0.5, 0.05)),
      kappa   = c(rbeta(n_per_chain, 20, 5),
                  rbeta(n_per_chain, 5, 20)),
      v       = rlnorm(n_draws, -1.5, 0.1),
      E_G     = rlnorm(n_draws, 6.0, 0.1),
      E0      = rlnorm(n_draws, 0.0, 0.1),
      L0      = rlnorm(n_draws, -2.0, 0.1),
      sigma_L = abs(rnorm(n_draws, 0.015, 0.003))
    )
  } else {
    pars <- data.frame(
      p_Am    = rlnorm(n_draws, 1.5, 0.1),
      p_M     = rlnorm(n_draws, -0.7, 0.1),
      kappa   = rbeta(n_draws, 8, 3),
      v       = rlnorm(n_draws, -1.5, 0.1),
      E_G     = rlnorm(n_draws, 6.0, 0.1),
      E0      = rlnorm(n_draws, 0.0, 0.1),
      L0      = rlnorm(n_draws, -2.0, 0.1),
      sigma_L = abs(rnorm(n_draws, 0.015, 0.003))
    )
  }

  n_obs <- 5
  for (i in seq_len(n_obs)) {
    pars[[paste0("L_hat[", i, "]")]] <- 0.1 + 0.04 * i + rnorm(n_draws, 0, 0.005)
    pars[[paste0("L_rep[", i, "]")]] <- 0.1 + 0.04 * i + rnorm(n_draws, 0, 0.02)
    pars[[paste0("log_lik[", i, "]")]] <- rnorm(n_draws, -1, 0.5)
  }
  pars$lp__ <- rnorm(n_draws, -50, 5)

  pars$.chain     <- rep(seq_len(n_chains), each = n_per_chain)
  pars$.iteration <- rep(seq_len(n_per_chain), n_chains)
  pars$.draw      <- seq_len(n_draws)

  draws_obj <- posterior::as_draws_df(pars)

  mock_cmdstan <- list(
    draws = function(format = NULL, ...) {
      if (!is.null(format) && format == "draws_matrix") {
        return(posterior::as_draws_matrix(draws_obj))
      }
      draws_obj
    },
    diagnostic_summary = function(quiet = FALSE) {
      list(
        num_divergent     = rep(as.integer(n_divergent), n_chains),
        num_max_treedepth = rep(as.integer(n_max_treedepth), n_chains),
        ebfmi             = if (low_ebfmi) rep(0.15, n_chains)
                            else rep(1.1, n_chains)
      )
    }
  )

  stan_data <- list(
    N_obs  = n_obs,
    t_obs  = seq(0, 28, length.out = n_obs),
    L_obs  = 0.1 + 0.04 * seq_len(n_obs),
    f_food = 1.0
  )

  mock_data <- structure(list(
    growth    = data.frame(id = 1, time = stan_data$t_obs,
                           length = stan_data$L_obs),
    n_ind     = 1L, ids = "1", endpoints = "growth",
    f_food    = 1.0, concentration = NULL
  ), class = "bdeb_data")

  mock_model <- structure(list(
    data = mock_data, type = "individual",
    priors = prior_default("individual"),
    observation = list(growth = obs_normal(), reproduction = obs_negbinom()),
    temperature = NULL, stan_model_name = "bdeb_individual_growth",
    stan_data = stan_data
  ), class = "bdeb_model")

  structure(list(
    fit = mock_cmdstan, model = mock_model, stan_model = NULL,
    chains = n_chains, iter_warmup = 500L,
    iter_sampling = n_per_chain, adapt_delta = 0.8
  ), class = "bdeb_fit")
}
