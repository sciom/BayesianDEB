## =============================================================
## SBC validation for BayesianDEB DEBtox model
## =============================================================
##
## Simulation-based calibration (Talts et al. 2018) for the
## DEBtox (TKTD) model with assimilation stress.  Each replication:
##   1. Draws parameters from the prior
##   2. Forward-simulates a 4-state DEBtox (E, V, R, Dw) via LSODA
##   3. Generates Gaussian growth observations for 4 conc. groups
##   4. Fits the DEBtox model
##   5. Computes parameter ranks
##
## Output:
##   sbc_debtox_results.rds  — archived results
##   (console)               — chi-squared tests & coverage
##
## Runtime: ~6–12 hours for 100 replications on 16-core machine
## =============================================================

library(BayesianDEB)
library(posterior)

set.seed(2026)

## --- Configuration ---------------------------------------------------

N_REP      <- 100
L_DRAWS    <- 199
N_CHAINS   <- 4
N_WARMUP   <- 500
N_SAMPLING <- 500
ADAPT_DELTA <- 0.95
T_MAX      <- 84
T_OBS      <- seq(0, T_MAX, by = 7)  # 13 weekly observations
CONC_LEVELS <- c(0, 10, 50, 200)     # 4 concentration groups

## --- Priors (tighter than defaults for SBC feasibility) --------------

sbc_priors <- list(
  p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
  p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
  kappa   = prior_beta(a = 3, b = 2),
  v       = prior_lognormal(mu = -1.5, sigma = 0.5),
  E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
  E0      = prior_lognormal(mu = 0.0, sigma = 0.5),
  L0      = prior_lognormal(mu = -2.0, sigma = 0.5),
  sigma_L = prior_halfnormal(sigma = 0.05),
  k_d     = prior_lognormal(mu = -1.0, sigma = 0.5),
  z_w     = prior_lognormal(mu = 3.0, sigma = 0.5),
  b_w     = prior_lognormal(mu = -4.0, sigma = 0.5)
)

## --- Draw from prior -------------------------------------------------

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
    k_d     = rlnorm(1, -1.0, 0.5),
    z_w     = rlnorm(1, 3.0, 0.5),
    b_w     = rlnorm(1, -4.0, 0.5)
  )
}

## --- Simulate DEBtox data for all concentrations ---------------------

simulate_debtox_data <- function(theta, t_obs, conc_levels) {
  all_data <- list()
  for (i in seq_along(conc_levels)) {
    C_w <- conc_levels[i]
    sim <- tryCatch(
      debtox_simulate(
        t_max = max(t_obs), p_Am = theta$p_Am, p_M = theta$p_M,
        kappa = theta$kappa, v = theta$v, E_G = theta$E_G,
        E0 = theta$E0, L0 = theta$L0, k_d = theta$k_d,
        z_w = theta$z_w, b_w = theta$b_w, C_w = C_w,
        f = 1.0, dt = 0.5),
      error = function(e) NULL
    )
    if (is.null(sim)) return(NULL)

    L_true <- approx(sim$time, sim$L, xout = t_obs)$y
    if (any(!is.finite(L_true)) || any(L_true <= 0)) return(NULL)

    L_obs <- rnorm(length(t_obs), L_true, theta$sigma_L)
    if (any(L_obs <= 0)) return(NULL)

    all_data[[i]] <- data.frame(
      id = i, time = t_obs, length = L_obs,
      concentration = C_w
    )
  }
  do.call(rbind, all_data)
}

## --- SBC parameters to track ----------------------------------------

sbc_pars <- c("p_Am", "p_M", "kappa", "v", "E_G",
              "sigma_L", "k_d", "z_w", "b_w")

## --- SBC loop --------------------------------------------------------

cat("=== SBC for DEBtox model ===\n")
cat(sprintf("Replications: %d, Concentrations: %s\n",
            N_REP, paste(CONC_LEVELS, collapse = ", ")))

results <- list()
n_valid <- 0
n_degenerate <- 0
n_divergent <- 0

for (r in seq_len(N_REP)) {
  cat(sprintf("\n--- Replication %d/%d ---\n", r, N_REP))

  ## 1. Draw from prior
  theta_star <- draw_prior()

  ## 2. Forward-simulate
  df <- simulate_debtox_data(theta_star, T_OBS, CONC_LEVELS)
  if (is.null(df)) {
    cat("  SKIP: degenerate simulation\n")
    n_degenerate <- n_degenerate + 1
    next
  }

  ## 3. Fit model
  conc_map <- setNames(CONC_LEVELS, as.character(CONC_LEVELS))
  dat <- tryCatch(
    bdeb_data(growth = df, concentration = conc_map, f_food = 1.0),
    error = function(e) NULL
  )
  if (is.null(dat)) {
    cat("  SKIP: data preparation failed\n")
    n_degenerate <- n_degenerate + 1
    next
  }

  mod <- tryCatch(
    bdeb_tox(dat, stress = "assimilation", priors = sbc_priors),
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
             threads_per_chain = 2,
             refresh = 0, seed = r),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    cat("  SKIP: fitting failed\n")
    n_divergent <- n_divergent + 1
    next
  }

  ## 4. Check convergence
  diag <- fit$fit$diagnostic_summary(quiet = TRUE)
  n_div <- sum(diag$num_divergent)

  draws <- as_draws_matrix(fit$fit$draws())
  avail_pars <- intersect(sbc_pars, colnames(draws))
  rhat_vals <- summarise_draws(draws[, avail_pars],
    rhat = posterior::rhat)$rhat

  if (n_div > 0 || any(rhat_vals > 1.05, na.rm = TRUE)) {
    cat(sprintf("  SKIP: divergences=%d, max_Rhat=%.3f\n",
                n_div, max(rhat_vals, na.rm = TRUE)))
    n_divergent <- n_divergent + 1
    next
  }

  ## 5. Compute ranks
  draws_mat <- as.data.frame(draws)
  idx <- sort(sample.int(nrow(draws_mat),
                         min(L_DRAWS, nrow(draws_mat))))
  draws_thin <- draws_mat[idx, avail_pars, drop = FALSE]

  true_vals <- c(
    p_Am = theta_star$p_Am, p_M = theta_star$p_M,
    kappa = theta_star$kappa, v = theta_star$v,
    E_G = theta_star$E_G, sigma_L = theta_star$sigma_L,
    k_d = theta_star$k_d, z_w = theta_star$z_w,
    b_w = theta_star$b_w
  )

  ranks <- sapply(avail_pars, function(p) {
    sum(draws_thin[[p]] < true_vals[p])
  })

  ## Also track derived EC50
  if ("z_w" %in% colnames(draws_mat) && "b_w" %in% colnames(draws_mat)) {
    ec50_post <- draws_mat[idx, "z_w"] + 0.5 / draws_mat[idx, "b_w"]
    ec50_true <- theta_star$z_w + 0.5 / theta_star$b_w
    ranks["EC50"] <- sum(ec50_post < ec50_true)
  }

  n_valid <- n_valid + 1
  results[[n_valid]] <- list(
    rep = r, ranks = ranks, true = true_vals,
    rhat = rhat_vals, n_div = n_div
  )

  cat(sprintf("  VALID (%d total)\n", n_valid))
}

## --- Summary ---------------------------------------------------------

cat(sprintf("\n=== SBC Summary ===\n"))
cat(sprintf("Total: %d, Valid: %d, Degenerate: %d, Divergent: %d\n",
            N_REP, n_valid, n_degenerate, n_divergent))

if (n_valid >= 20) {
  all_tracked <- names(results[[1]]$ranks)
  rank_mat <- do.call(rbind, lapply(results, function(x) x$ranks))

  n_bins <- min(20, floor(n_valid / 5))
  cat(sprintf("\nChi-squared tests (%d bins, %d df):\n",
              n_bins, n_bins - 1))

  for (p in all_tracked) {
    h <- hist(rank_mat[, p], breaks = seq(0, L_DRAWS,
              length.out = n_bins + 1), plot = FALSE)
    expected <- n_valid / n_bins
    chi2 <- sum((h$counts - expected)^2 / expected)
    pval <- 1 - pchisq(chi2, df = n_bins - 1)
    cat(sprintf("  %-10s: chi2 = %5.1f, p = %.3f\n", p, chi2, pval))
  }

  cat("\nCredible interval coverage:\n")
  cat(sprintf("  %-10s  50%%    90%%    95%%\n", "Parameter"))
  for (p in all_tracked) {
    cov50 <- mean(rank_mat[, p] >= L_DRAWS * 0.25 &
                  rank_mat[, p] <= L_DRAWS * 0.75)
    cov90 <- mean(rank_mat[, p] >= L_DRAWS * 0.05 &
                  rank_mat[, p] <= L_DRAWS * 0.95)
    cov95 <- mean(rank_mat[, p] >= L_DRAWS * 0.025 &
                  rank_mat[, p] <= L_DRAWS * 0.975)
    cat(sprintf("  %-10s  %.2f   %.2f   %.2f\n",
                p, cov50, cov90, cov95))
  }
}

## --- Save results ----------------------------------------------------

saveRDS(list(
  results = results, n_valid = n_valid, n_total = N_REP,
  n_degenerate = n_degenerate, n_divergent = n_divergent,
  sbc_pars = sbc_pars, L_draws = L_DRAWS,
  conc_levels = CONC_LEVELS,
  timestamp = Sys.time()
), "sbc_debtox_results.rds")

cat("\nResults saved to sbc_debtox_results.rds\n")
cat("=== DONE ===\n")
