# ****************************************************************
# SBC validation for BayesianDEB growth_repro model
# Uses E. fetida-like parameters (same as individual SBC) with
# added reproduction parameters.  SBC validates the inference
# pipeline, not species-specific biology.
# Output: sbc_growth_repro_results.rds (archived results)
# Runtime: cca 2–4 hours for 100 replications on 8-core machine
# ****************************************************************

library(BayesianDEB)
library(posterior)
library(deSolve)

set.seed(2026)

# ****************************************************************
# Configuration
# ****************************************************************

N_REP      <- 100   # total replications
L_DRAWS    <- 199   # thinned posterior draws per replication
N_CHAINS   <- 4
N_WARMUP   <- 500
N_SAMPLING <- 500
ADAPT_DELTA <- 0.9
T_MAX      <- 84
N_GROWTH   <- 13    # weekly growth observations

sbc_priors <- list(
  p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
  p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
  kappa   = prior_beta(a = 3, b = 2),
  v       = prior_lognormal(mu = -1.5, sigma = 0.5),
  E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
  E0      = prior_lognormal(mu = 0.0, sigma = 0.5),
  L0      = prior_lognormal(mu = -2.0, sigma = 0.5),
  sigma_L = prior_halfnormal(sigma = 0.05),
  k_J     = prior_lognormal(mu = -5.0, sigma = 0.5),
  E_Hp    = prior_lognormal(mu = -4.0, sigma = 0.5),
  k_R     = prior_lognormal(mu = -2.0, sigma = 0.5),
  phi_R   = prior_lognormal(mu = 2.0, sigma = 0.5)
)

# ****************************************************************
# ODE system for 3-state DEB (E, V, R)
# ****************************************************************

deb_3state_ode <- function(t, state, pars) {
  E <- state[1]; V <- state[2]; R <- state[3]
  if (V <= 0) return(list(c(0, 0, 0)))
  L <- V^(1/3)

  p_A <- pars$f * pars$p_Am * L^2
  p_C <- E * pars$v * L / (E + pars$E_G * V)
  p_M <- pars$p_M * V

  dE <- p_A - p_C
  dV <- max((pars$kappa * p_C - p_M) / pars$E_G, 0)
  dR <- max((1 - pars$kappa) * p_C - pars$k_J * pars$E_Hp, 0)

  list(c(dE, dV, dR))
}

# ****************************************************************
# Draw from prior
# ****************************************************************

draw_prior <- function() {
  list(
    p_Am    = rlnorm(1, 1.5, 0.5),
    p_M     = rlnorm(1, -1.0, 0.5),
    kappa   = rbeta(1, 3, 2),
    v       = rlnorm(1, -1.5, 0.5),
    E_G     = rlnorm(1, 6.0, 0.5),
    E0      = rlnorm(1, 0.0, 0.5),
    L0      = rlnorm(1, -2.0, 0.5),
    sigma_L = abs(rnorm(1, 0, 0.05)),
    k_J     = rlnorm(1, -5.0, 0.5),
    E_Hp    = rlnorm(1, -4.0, 0.5),
    k_R     = rlnorm(1, -2.0, 0.5),
    phi_R   = rlnorm(1, 2.0, 0.5),
    f       = 1.0
  )
}

# ****************************************************************
# Forward-simulate
# ****************************************************************

simulate_data <- function(theta, t_growth, t_repro_start, t_repro_end) {
  V0 <- theta$L0^3
  state0 <- c(E = theta$E0, V = V0, R = 0)

  # All unique times needed (sorted, no duplicates, no exact zeros)
  t_all <- sort(unique(c(1e-3, t_growth, t_repro_start, t_repro_end)))
  t_all <- t_all[t_all > 0]

  sol <- tryCatch(
    lsoda(state0, t_all, deb_3state_ode, theta,
           rtol = 1e-6, atol = 1e-6, maxsteps = 10000),
    error = function(e) NULL
  )
  if (is.null(sol) || nrow(sol) < length(t_all)) return(NULL)

  # Growth observations
  idx_growth <- match(t_growth, sol[, "time"])
  if (any(is.na(idx_growth))) return(NULL)
  L_true <- sol[idx_growth, "V"]^(1/3)
  if (any(!is.finite(L_true)) || any(L_true <= 0)) return(NULL)

  L_obs <- rnorm(length(t_growth), L_true, theta$sigma_L)
  if (any(L_obs <= 0)) return(NULL)

  # Reproduction: accumulated buffer over interval
  idx_start <- match(t_repro_start, sol[, "time"])
  idx_end   <- match(t_repro_end, sol[, "time"])
  if (any(is.na(c(idx_start, idx_end)))) return(NULL)

  delta_R <- sol[idx_end, "R"] - sol[idx_start, "R"]
  if (any(!is.finite(delta_R)) || any(delta_R < 0)) return(NULL)

  # Generate offspring counts
  mu_R <- theta$k_R * delta_R
  if (any(!is.finite(mu_R)) || any(mu_R <= 0)) return(NULL)

  R_counts <- rnbinom(length(mu_R), size = theta$phi_R, mu = mu_R)

  growth_df <- data.frame(id = 1, time = t_growth, length = L_obs)
  repro_df  <- data.frame(id = 1, t_start = t_repro_start,
                           t_end = t_repro_end, count = R_counts)

  list(growth = growth_df, repro = repro_df)
}

# ****************************************************************
# SBC parameters to track
# ****************************************************************

sbc_pars <- c("p_Am", "p_M", "kappa", "v", "E_G",
              "sigma_L", "k_R", "phi_R")

# --- SBC loop --------------------------------------------------------

cat("*** SBC for growth_repro model ***\n")
cat(sprintf("Replications: %d, Chains: %d\n", N_REP, N_CHAINS))

t_growth <- seq(0, T_MAX, length.out = N_GROWTH)
t_growth[1] <- 0.001  # avoid t=0
# Two reproduction intervals: 0-42 and 42-84 days
t_repro_start <- c(0.001, 42)
t_repro_end   <- c(42, 84)

results <- list()
n_valid <- 0
n_degenerate <- 0
n_divergent <- 0

for (r in seq_len(N_REP)) {
  cat(sprintf("\n--- Replication %d/%d ---\n", r, N_REP))

  # 1. Draw from prior
  theta_star <- draw_prior()

  # 2. Forward-simulate
  sim <- simulate_data(theta_star, t_growth, t_repro_start, t_repro_end)
  if (is.null(sim)) {
    cat("  SKIP: degenerate simulation\n")
    n_degenerate <- n_degenerate + 1
    next
  }

  # 3. Fit model
  dat <- tryCatch(
    bdeb_data(growth = sim$growth, reproduction = sim$repro,
              f_food = 1.0),
    error = function(e) NULL
  )
  if (is.null(dat)) {
    cat("  SKIP: data preparation failed\n")
    n_degenerate <- n_degenerate + 1
    next
  }

  mod <- tryCatch(
    bdeb_model(dat, type = "growth_repro",
      priors = sbc_priors,
      observation = list(growth = obs_normal(),
                         reproduction = obs_negbinom())),
    error = function(e) NULL
  )
  if (is.null(mod)) {
    cat("  SKIP: model specification failed\n")
    n_degenerate <- n_degenerate + 1
    next
  }

  fit <- tryCatch(
    bdeb_fit(mod, chains = N_CHAINS,
             iter_warmup = N_WARMUP,
             iter_sampling = N_SAMPLING,
             adapt_delta = ADAPT_DELTA,
             refresh = 0, seed = r),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    cat("  SKIP: fitting failed\n")
    n_divergent <- n_divergent + 1
    next
  }

  # 4. Check convergence
  diag <- fit$fit$diagnostic_summary(quiet = TRUE)
  n_div <- sum(diag$num_divergent)

  draws <- as_draws_matrix(fit$fit$draws())
  avail_pars <- intersect(sbc_pars, colnames(draws))
  if (length(avail_pars) == 0) {
    cat("  SKIP: no parameters found\n")
    n_divergent <- n_divergent + 1
    next
  }

  rhat_vals <- summarise_draws(draws[, avail_pars],
    rhat = posterior::rhat)$rhat

  if (n_div > 0 || any(rhat_vals > 1.05, na.rm = TRUE)) {
    cat(sprintf("  SKIP: divergences=%d, max_Rhat=%.3f\n",
                n_div, max(rhat_vals, na.rm = TRUE)))
    n_divergent <- n_divergent + 1
    next
  }

  # 5. Compute ranks
  draws_mat <- as.data.frame(draws)
  idx <- sort(sample.int(nrow(draws_mat),
                         min(L_DRAWS, nrow(draws_mat))))
  draws_thin <- draws_mat[idx, avail_pars, drop = FALSE]

  true_vals <- c(
    p_Am = theta_star$p_Am, p_M = theta_star$p_M,
    kappa = theta_star$kappa, v = theta_star$v,
    E_G = theta_star$E_G, sigma_L = theta_star$sigma_L,
    k_R = theta_star$k_R, phi_R = theta_star$phi_R
  )

  ranks <- sapply(avail_pars, function(p) {
    sum(draws_thin[[p]] < true_vals[p])
  })

  n_valid <- n_valid + 1
  results[[n_valid]] <- list(
    rep = r, ranks = ranks, true = true_vals[avail_pars],
    rhat = rhat_vals, n_div = n_div
  )

  cat(sprintf("  VALID (%d total)\n", n_valid))
}

# ****************************************************************
# Summary
# ****************************************************************

cat(sprintf("\n*** SBC Summary ***\n"))
cat(sprintf("Total: %d, Valid: %d, Degenerate: %d, Divergent: %d\n",
            N_REP, n_valid, n_degenerate, n_divergent))

if (n_valid >= 10) {
  rank_mat <- do.call(rbind, lapply(results, function(x) x$ranks))

  n_bins <- min(20, max(5, floor(n_valid / 5)))
  cat(sprintf("\nChi-squared tests (%d bins, %d df):\n",
              n_bins, n_bins - 1))

  for (p in colnames(rank_mat)) {
    h <- hist(rank_mat[, p], breaks = seq(0, L_DRAWS,
              length.out = n_bins + 1), plot = FALSE)
    expected <- n_valid / n_bins
    chi2 <- sum((h$counts - expected)^2 / expected)
    pval <- 1 - pchisq(chi2, df = n_bins - 1)
    cat(sprintf("  %-10s: chi2 = %5.1f, p = %.3f\n", p, chi2, pval))
  }

  cat("\nCredible interval coverage:\n")
  cat(sprintf("  %-10s  50%%    90%%    95%%\n", "Parameter"))
  for (p in colnames(rank_mat)) {
    cov50 <- mean(rank_mat[, p] >= L_DRAWS * 0.25 &
                  rank_mat[, p] <= L_DRAWS * 0.75)
    cov90 <- mean(rank_mat[, p] >= L_DRAWS * 0.05 &
                  rank_mat[, p] <= L_DRAWS * 0.95)
    cov95 <- mean(rank_mat[, p] >= L_DRAWS * 0.025 &
                  rank_mat[, p] <= L_DRAWS * 0.975)
    cat(sprintf("  %-10s  %.2f   %.2f   %.2f\n",
                p, cov50, cov90, cov95))
  }
} else {
  cat("Too few valid replications for chi-squared testing.\n")
}

# --- Save results ----------------------------------------------------

saveRDS(list(
  results = results, n_valid = n_valid, n_total = N_REP,
  n_degenerate = n_degenerate, n_divergent = n_divergent,
  sbc_pars = sbc_pars, L_draws = L_DRAWS,
  timestamp = Sys.time()
), "sbc_growth_repro_results.rds")

cat("\nResults saved to sbc_growth_repro_results.rds\n")
cat("*** DONE ***\n")
