## =============================================================
## Simulation-Based Calibration (SBC) for BayesianDEB
## Reference: Talts et al. (2020), arXiv:1804.06788
## =============================================================
##
## Validates the individual growth model by checking that posterior
## ranks of prior-drawn parameters are uniformly distributed.
##
## Sequential loop with parallel_chains=4 inside each fit.
## Run time: ~45 min on a multi-core machine.
##
## Output:
##   fig_sbc_ranks.pdf   - rank histogram figure for manuscript
##   sbc_results.rds     - full results object
##   (console)           - LaTeX table rows for bayesiandeb.tex
## =============================================================

library(BayesianDEB)
library(posterior)
library(ggplot2)

## --- Configuration ---------------------------------------------------

N_SBC         <- 500L       # SBC replications
L_THIN        <- 199L       # posterior draws per replication (thinned)
CHAINS        <- 4L
ITER_WARMUP   <- 500L
ITER_SAMPLING <- 500L       # per chain; total raw = 2000
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 8L         # keep low to avoid stalling on hard params
SEED_BASE     <- 20260407L
DT_SIM        <- 0.05       # Euler step for data generation (days)
T_OBS         <- seq(0, 84, by = 7)  # 13 weekly time points
F_FOOD        <- 1.0
N_BINS        <- 20L        # bins for chi-squared test
FIT_TIMEOUT   <- 120        # seconds: skip fit if longer than this

## --- Priors (shared between simulation and fitting) ------------------

sbc_priors <- list(
  p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
  p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
  kappa   = prior_beta(a = 3, b = 2),
  v       = prior_lognormal(mu = -1.5, sigma = 0.5),
  E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
  E0      = prior_lognormal(mu = 0.0, sigma = 0.5),
  L0      = prior_lognormal(mu = -2.0, sigma = 0.5),
  sigma_L = prior_halfnormal(sigma = 0.05)
)

par_names <- c("p_Am", "p_M", "kappa", "v", "E_G", "E0", "L0", "sigma_L")
n_par     <- length(par_names)

## --- Prior sampling function -----------------------------------------

sample_prior <- function() {
  list(
    p_Am    = rlnorm(1, meanlog = 1.5,  sdlog = 0.5),
    p_M     = rlnorm(1, meanlog = -1.0, sdlog = 0.5),
    kappa   = rbeta(1, shape1 = 3, shape2 = 2),
    v       = rlnorm(1, meanlog = -1.5, sdlog = 0.5),
    E_G     = rlnorm(1, meanlog = 6.0,  sdlog = 0.5),
    E0      = rlnorm(1, meanlog = 0.0,  sdlog = 0.5),
    L0      = rlnorm(1, meanlog = -2.0, sdlog = 0.5),
    sigma_L = abs(rnorm(1, mean = 0, sd = 0.05))
  )
}

## --- Data simulation -------------------------------------------------

simulate_sbc_data <- function(theta, t_obs, f = 1.0, dt = 0.05) {
  sim <- deb_simulate(
    t_max = max(t_obs),
    p_Am  = theta$p_Am, p_M = theta$p_M, kappa = theta$kappa,
    v = theta$v, E_G = theta$E_G, E0 = theta$E0, L0 = theta$L0,
    f = f, dt = dt
  )
  L_true <- approx(sim$time, sim$L, xout = t_obs)$y
  if (any(is.na(L_true)) || any(L_true <= 0)) return(NULL)
  if (max(L_true) < 1e-6) return(NULL)
  L_obs <- rnorm(length(t_obs), mean = L_true, sd = theta$sigma_L)
  if (any(L_obs <= 0)) return(NULL)
  data.frame(id = 1, time = t_obs, length = L_obs)
}

## --- Pre-compile Stan model (one-time cost) --------------------------

cat("Pre-compiling Stan model...\n")
dummy_df <- data.frame(
  id = 1, time = T_OBS,
  length = seq(0.1, 1.0, length.out = length(T_OBS))
)
dummy_dat <- bdeb_data(growth = dummy_df, f_food = F_FOOD)
dummy_mod <- bdeb_model(dummy_dat, type = "individual", priors = sbc_priors)
dummy_fit <- bdeb_fit(dummy_mod, chains = 1, iter_warmup = 5,
                      iter_sampling = 5, refresh = 0, seed = 1)
cat("Compilation complete.\n\n")

## --- SBC main loop (sequential) --------------------------------------

ranks      <- matrix(NA_integer_, nrow = N_SBC, ncol = n_par,
                     dimnames = list(NULL, par_names))
converged  <- rep(FALSE, N_SBC)
theta_true <- matrix(NA_real_, nrow = N_SBC, ncol = n_par,
                     dimnames = list(NULL, par_names))
cov_levels <- c("50", "80", "90", "95")
coverage   <- array(FALSE, dim = c(N_SBC, n_par, length(cov_levels)),
                    dimnames = list(NULL, par_names, cov_levels))

set.seed(SEED_BASE)
seeds <- sample.int(1e7, N_SBC)

cat(sprintf("=== Starting SBC: %d replications ===\n", N_SBC))
t_start <- Sys.time()
n_skip <- 0L

for (i in seq_len(N_SBC)) {
  set.seed(seeds[i])

  ## 1. Draw from prior
  theta <- sample_prior()
  for (p in par_names) theta_true[i, p] <- theta[[p]]

  ## 2. Simulate data
  df <- simulate_sbc_data(theta, T_OBS, f = F_FOOD, dt = DT_SIM)
  if (is.null(df)) { n_skip <- n_skip + 1L; next }

  ## 3. Fit and compute ranks
  fit_ok <- tryCatch({
    dat <- bdeb_data(growth = df, f_food = F_FOOD)
    mod <- bdeb_model(dat, type = "individual", priors = sbc_priors)

    t_fit <- system.time({
      fit <- bdeb_fit(mod, chains = CHAINS, iter_warmup = ITER_WARMUP,
                      iter_sampling = ITER_SAMPLING, adapt_delta = ADAPT_DELTA,
                      max_treedepth = MAX_TREEDEPTH, parallel_chains = CHAINS,
                      refresh = 0, seed = seeds[i])
    })[["elapsed"]]

    # Skip if fit took too long (ODE solver struggles)
    if (t_fit > FIT_TIMEOUT) return(FALSE)

    # Check convergence
    diag <- bdeb_diagnose(fit)
    if (diag$n_divergent > 0) return(FALSE)

    smry <- bdeb_summary(fit, pars = par_names)
    if (any(smry$rhat > 1.05, na.rm = TRUE)) return(FALSE)

    # Extract and thin posterior draws
    draws_mat <- as_draws_matrix(fit$fit$draws())
    n_total   <- nrow(draws_mat)
    idx       <- sort(sample.int(n_total, min(L_THIN, n_total)))

    # Compute ranks and coverage
    for (j in seq_along(par_names)) {
      p    <- par_names[j]
      post <- as.numeric(draws_mat[idx, p])
      ranks[i, p] <- sum(post < theta[[p]])

      qs <- quantile(post, probs = c(0.025, 0.05, 0.10, 0.25,
                                      0.75, 0.90, 0.95, 0.975))
      tv <- theta[[p]]
      coverage[i, j, "50"] <- tv >= qs["25%"]  && tv <= qs["75%"]
      coverage[i, j, "80"] <- tv >= qs["10%"]  && tv <= qs["90%"]
      coverage[i, j, "90"] <- tv >= qs["5%"]   && tv <= qs["95%"]
      coverage[i, j, "95"] <- tv >= qs["2.5%"] && tv <= qs["97.5%"]
    }
    TRUE
  }, error = function(e) FALSE)

  converged[i] <- isTRUE(fit_ok)

  if (i %% 25 == 0) {
    n_ok    <- sum(converged)
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    eta     <- elapsed / i * (N_SBC - i)
    cat(sprintf("[%3d/%d] valid: %d | skip: %d | elapsed: %.0f min | ETA: %.0f min\n",
                i, N_SBC, n_ok, n_skip, elapsed, eta))
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nSBC complete in %.1f minutes\n\n", elapsed_total))

## --- Analyse results -------------------------------------------------

valid   <- which(converged)
N_valid <- length(valid)
cat(sprintf("Valid replications: %d / %d (%.0f%%)\n",
            N_valid, N_SBC, 100 * N_valid / N_SBC))
cat(sprintf("Skipped (degenerate data): %d\n", n_skip))
cat(sprintf("Failed (divergences/Rhat/timeout): %d\n",
            N_SBC - N_valid - n_skip))

ranks_valid    <- ranks[valid, , drop = FALSE]
coverage_valid <- coverage[valid, , , drop = FALSE]

## --- Chi-squared uniformity test -------------------------------------

bin_breaks     <- seq(0, L_THIN, length.out = N_BINS + 1)
expected_count <- N_valid / N_BINS

sbc_table <- data.frame(
  parameter = par_names,
  chi_sq    = NA_real_,
  p_value   = NA_real_,
  cov_50    = NA_real_,
  cov_80    = NA_real_,
  cov_90    = NA_real_,
  cov_95    = NA_real_,
  stringsAsFactors = FALSE
)

for (j in seq_along(par_names)) {
  p <- par_names[j]
  r <- ranks_valid[, p]

  h   <- hist(r, breaks = bin_breaks, plot = FALSE)
  chi <- chisq.test(h$counts)
  sbc_table$chi_sq[j]  <- unname(chi$statistic)
  sbc_table$p_value[j] <- chi$p.value

  sbc_table$cov_50[j] <- mean(coverage_valid[, j, "50"])
  sbc_table$cov_80[j] <- mean(coverage_valid[, j, "80"])
  sbc_table$cov_90[j] <- mean(coverage_valid[, j, "90"])
  sbc_table$cov_95[j] <- mean(coverage_valid[, j, "95"])
}

cat("\n=== SBC Summary ===\n")
print(sbc_table, digits = 3, row.names = FALSE)

## --- LaTeX table output ----------------------------------------------

cat("\n=== LaTeX table rows (paste into bayesiandeb.tex, tab:sbc) ===\n\n")

tex_names <- c("$\\pAm$", "$\\pM$", "$\\kappa$", "$v$",
               "$\\EG$", "$E_0$", "$L_0$", "$\\sigma_L$")

show_idx <- c(1, 2, 3, 4, 5, 8)  # p_Am, p_M, kappa, v, E_G, sigma_L
for (j in show_idx) {
  cat(sprintf("%s & %.1f & %.2f & %.2f & %.2f & %.2f \\\\\n",
              tex_names[j],
              sbc_table$chi_sq[j], sbc_table$p_value[j],
              sbc_table$cov_50[j], sbc_table$cov_90[j],
              sbc_table$cov_95[j]))
}
cat(sprintf("\n%% N_valid = %d / %d\n", N_valid, N_SBC))

## --- Save full results (before plotting, to avoid data loss) ---------

saveRDS(list(
  ranks       = ranks_valid,
  coverage    = coverage_valid,
  theta_true  = theta_true[valid, ],
  sbc_table   = sbc_table,
  converged   = converged,
  config      = list(N_SBC = N_SBC, N_valid = N_valid,
                     L_THIN = L_THIN, CHAINS = CHAINS,
                     ITER_WARMUP = ITER_WARMUP,
                     ITER_SAMPLING = ITER_SAMPLING,
                     ADAPT_DELTA = ADAPT_DELTA,
                     seeds = seeds, priors = sbc_priors)
), "sbc_results.rds")
cat("Results saved to sbc_results.rds\n")

## --- Rank histogram figure -------------------------------------------

show_pars   <- c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L")
show_labels <- c(
  p_Am    = "p[Am]",
  p_M     = "p[M]",
  kappa   = "kappa",
  v       = "v",
  E_G     = "E[G]",
  sigma_L = "sigma[L]"
)

ranks_long <- do.call(rbind, lapply(show_pars, function(p) {
  data.frame(parameter = show_labels[p],
             rank = ranks_valid[, p],
             stringsAsFactors = FALSE)
}))
ranks_long$parameter <- factor(ranks_long$parameter,
                               levels = unname(show_labels))

p_sbc <- ggplot(ranks_long, aes(x = rank)) +
  geom_histogram(bins = N_BINS, fill = "steelblue", colour = "white",
                 linewidth = 0.3, alpha = 0.85) +
  geom_hline(yintercept = expected_count, linetype = "dashed",
             colour = "firebrick", linewidth = 0.5) +
  facet_wrap(~ parameter, ncol = 3, labeller = label_parsed,
             scales = "fixed") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 11)) +
  labs(x = "Rank", y = "Count")

ggsave("fig_sbc_ranks.pdf", p_sbc, width = 7, height = 4.5)
cat("\nRank histogram saved to fig_sbc_ranks.pdf\n")

cat("=== DONE ===\n")
