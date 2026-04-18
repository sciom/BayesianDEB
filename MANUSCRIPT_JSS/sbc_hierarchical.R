## =============================================================
## SBC for Hierarchical Growth Model (BayesianDEB)
## 50 replications, 5 individuals, 13 time points
## Optimised: 2 chains x 8 threads, 250 iter
## =============================================================

library(BayesianDEB)
library(posterior)
library(ggplot2)

## --- Configuration ---------------------------------------------------

N_SBC         <- 50L
L_THIN        <- 99L          # thinned draws (250 sampling / 2 chains ~ 500 total)
CHAINS        <- 2L           # fewer chains, more threads
THREADS       <- 8L           # threads per chain (reduce_sum)
ITER_WARMUP   <- 250L         # reduced from 500
ITER_SAMPLING <- 250L         # reduced from 500
ADAPT_DELTA   <- 0.95
MAX_TREEDEPTH <- 10L
SEED_BASE     <- 20260408L
DT_SIM        <- 0.05
N_IND         <- 5L           # reduced from 10
T_OBS         <- seq(0, 84, by = 7)  # 13 weekly
F_FOOD        <- 1.0
N_BINS        <- 10L
FIT_TIMEOUT   <- 180          # seconds (reduced from 300)

## --- Priors (must match bdeb_model specification) --------------------

sbc_priors <- list(
  mu_log_p_Am    = prior_normal(mu = 1.5, sigma = 0.5),
  sigma_log_p_Am = prior_exponential(rate = 2),
  p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
  kappa   = prior_beta(a = 3, b = 2),
  v       = prior_lognormal(mu = -1.5, sigma = 0.5),
  E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
  E0      = prior_lognormal(mu = 0.0, sigma = 0.5),
  L0      = prior_lognormal(mu = -2.0, sigma = 0.5),
  sigma_L = prior_halfnormal(sigma = 0.05)
)

## --- Prior sampling --------------------------------------------------

sample_hier_prior <- function(n_ind) {
  mu_log_pAm    <- rnorm(1, 1.5, 0.5)
  sigma_log_pAm <- rexp(1, 2)
  p_Am_ind <- exp(rnorm(n_ind, mu_log_pAm, sigma_log_pAm))

  list(
    mu_log_p_Am    = mu_log_pAm,
    sigma_log_p_Am = sigma_log_pAm,
    p_Am_ind = p_Am_ind,
    p_M     = rlnorm(1, -1.0, 0.5),
    kappa   = rbeta(1, 3, 2),
    v       = rlnorm(1, -1.5, 0.5),
    E_G     = rlnorm(1, 6.0, 0.5),
    E0      = rlnorm(1, 0.0, 0.5),
    L0      = rlnorm(1, -2.0, 0.5),
    sigma_L = abs(rnorm(1, 0, 0.05))
  )
}

## --- Data simulation -------------------------------------------------

simulate_hier_data <- function(theta, n_ind, t_obs, f = 1.0, dt = 0.05) {
  rows <- list()
  for (j in seq_len(n_ind)) {
    sim <- deb_simulate(
      t_max = max(t_obs), p_Am = theta$p_Am_ind[j],
      p_M = theta$p_M, kappa = theta$kappa, v = theta$v,
      E_G = theta$E_G, E0 = theta$E0, L0 = theta$L0,
      f = f, dt = dt
    )
    L_true <- approx(sim$time, sim$L, xout = t_obs)$y
    if (any(is.na(L_true)) || any(L_true <= 0)) return(NULL)
    L_obs <- rnorm(length(t_obs), L_true, theta$sigma_L)
    if (any(L_obs <= 0)) return(NULL)
    rows[[j]] <- data.frame(id = j, time = t_obs, length = L_obs)
  }
  do.call(rbind, rows)
}

## --- Pre-compile Stan model ------------------------------------------

cat("Pre-compiling hierarchical Stan model...\n")
dummy_df <- data.frame(
  id = rep(1:N_IND, each = length(T_OBS)),
  time = rep(T_OBS, N_IND),
  length = rep(seq(0.1, 0.5, length.out = length(T_OBS)), N_IND)
)
dummy_dat <- bdeb_data(growth = dummy_df, f_food = F_FOOD)
dummy_mod <- bdeb_model(dummy_dat, type = "hierarchical", priors = sbc_priors)
dummy_fit <- bdeb_fit(dummy_mod, chains = 1, iter_warmup = 5,
                      iter_sampling = 5, threads_per_chain = THREADS,
                      refresh = 0, seed = 1)
cat("Compilation complete.\n\n")

## --- SBC main loop ---------------------------------------------------

par_names <- c("mu_log_p_Am", "sigma_log_p_Am", "p_M", "kappa", "sigma_L")
n_par     <- length(par_names)

ranks      <- matrix(NA_integer_, nrow = N_SBC, ncol = n_par,
                     dimnames = list(NULL, par_names))
converged  <- rep(FALSE, N_SBC)
theta_true <- matrix(NA_real_, nrow = N_SBC, ncol = n_par,
                     dimnames = list(NULL, par_names))

cov_levels <- c("50", "90", "95")
coverage   <- array(FALSE, dim = c(N_SBC, n_par, length(cov_levels)),
                    dimnames = list(NULL, par_names, cov_levels))

set.seed(SEED_BASE)
seeds <- sample.int(1e7, N_SBC)

cat(sprintf("=== Starting hierarchical SBC: %d replications, %d individuals ===\n",
            N_SBC, N_IND))
t_start <- Sys.time()
n_skip <- 0L

for (i in seq_len(N_SBC)) {
  set.seed(seeds[i])

  theta <- sample_hier_prior(N_IND)
  for (p in par_names) theta_true[i, p] <- theta[[p]]

  df <- simulate_hier_data(theta, N_IND, T_OBS, f = F_FOOD, dt = DT_SIM)
  if (is.null(df)) { n_skip <- n_skip + 1L; next }

  fit_ok <- tryCatch({
    dat <- bdeb_data(growth = df, f_food = F_FOOD)
    mod <- bdeb_model(dat, type = "hierarchical", priors = sbc_priors)

    t_fit <- system.time({
      fit <- bdeb_fit(mod, chains = CHAINS, iter_warmup = ITER_WARMUP,
                      iter_sampling = ITER_SAMPLING, adapt_delta = ADAPT_DELTA,
                      max_treedepth = MAX_TREEDEPTH, parallel_chains = CHAINS,
                      threads_per_chain = THREADS,
                      refresh = 0, seed = seeds[i])
    })[["elapsed"]]

    if (t_fit > FIT_TIMEOUT) return(FALSE)

    diag <- bdeb_diagnose(fit)
    if (diag$n_divergent > 0) return(FALSE)

    smry <- bdeb_summary(fit, pars = par_names)
    if (any(smry$rhat > 1.05, na.rm = TRUE)) return(FALSE)

    draws_mat <- as.data.frame(as_draws_matrix(fit$fit$draws()))
    n_total   <- nrow(draws_mat)
    idx       <- sort(sample.int(n_total, min(L_THIN, n_total)))

    for (j in seq_along(par_names)) {
      p    <- par_names[j]
      post <- as.numeric(draws_mat[idx, p])
      ranks[i, p] <- sum(post < theta[[p]])

      qs <- quantile(post, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
      tv <- theta[[p]]
      coverage[i, j, "50"] <- tv >= qs["25%"]  && tv <= qs["75%"]
      coverage[i, j, "90"] <- tv >= qs["5%"]   && tv <= qs["95%"]
      coverage[i, j, "95"] <- tv >= qs["2.5%"] && tv <= qs["97.5%"]
    }
    TRUE
  }, error = function(e) FALSE)

  converged[i] <- isTRUE(fit_ok)

  if (i %% 10 == 0) {
    n_ok    <- sum(converged)
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    eta     <- elapsed / i * (N_SBC - i)
    cat(sprintf("[%3d/%d] valid: %d | skip: %d | %.0f min | ETA: %.0f min\n",
                i, N_SBC, n_ok, n_skip, elapsed, eta))
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nSBC complete in %.1f minutes\n", elapsed_total))

## --- Results ---------------------------------------------------------

valid   <- which(converged)
N_valid <- length(valid)
cat(sprintf("Valid: %d / %d (%.0f%%)\n", N_valid, N_SBC, 100 * N_valid / N_SBC))

ranks_valid    <- ranks[valid, , drop = FALSE]
coverage_valid <- coverage[valid, , , drop = FALSE]

bin_breaks     <- seq(0, L_THIN, length.out = N_BINS + 1)
expected_count <- N_valid / N_BINS

sbc_table <- data.frame(parameter = par_names, chi_sq = NA_real_,
  p_value = NA_real_, cov_50 = NA_real_, cov_90 = NA_real_,
  cov_95 = NA_real_, stringsAsFactors = FALSE)

for (j in seq_along(par_names)) {
  r <- ranks_valid[, par_names[j]]
  h <- hist(r, breaks = bin_breaks, plot = FALSE)
  chi <- chisq.test(h$counts)
  sbc_table$chi_sq[j]  <- unname(chi$statistic)
  sbc_table$p_value[j] <- chi$p.value
  sbc_table$cov_50[j]  <- mean(coverage_valid[, j, "50"])
  sbc_table$cov_90[j]  <- mean(coverage_valid[, j, "90"])
  sbc_table$cov_95[j]  <- mean(coverage_valid[, j, "95"])
}

cat("\n=== SBC Summary ===\n")
print(sbc_table, digits = 3, row.names = FALSE)

## --- LaTeX output ----------------------------------------------------

cat("\n=== LaTeX rows ===\n")
tex_names <- c("$\\mu_{\\log \\pAm}$", "$\\sigma_{\\log \\pAm}$",
               "$\\pM$", "$\\kappa$", "$\\sigma_L$")
for (j in seq_along(par_names)) {
  cat(sprintf("%s & %.1f & %.2f & %.2f & %.2f & %.2f \\\\\n",
    tex_names[j], sbc_table$chi_sq[j], sbc_table$p_value[j],
    sbc_table$cov_50[j], sbc_table$cov_90[j], sbc_table$cov_95[j]))
}
cat(sprintf("\n%% N_valid = %d / %d\n", N_valid, N_SBC))

## --- Rank histogram --------------------------------------------------

show_labels <- c(
  mu_log_p_Am = "mu[log~p[Am]]",
  sigma_log_p_Am = "sigma[log~p[Am]]",
  p_M = "p[M]", kappa = "kappa", sigma_L = "sigma[L]"
)

ranks_long <- do.call(rbind, lapply(par_names, function(p) {
  data.frame(parameter = show_labels[p], rank = ranks_valid[, p],
             stringsAsFactors = FALSE)
}))
ranks_long$parameter <- factor(ranks_long$parameter, levels = unname(show_labels))

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

ggsave("fig_sbc_hierarchical.pdf", p_sbc, width = 7, height = 4.5)
cat("\nRank histogram saved to fig_sbc_hierarchical.pdf\n")

saveRDS(list(ranks = ranks_valid, coverage = coverage_valid,
  sbc_table = sbc_table, converged = converged,
  config = list(N_SBC = N_SBC, N_valid = N_valid, N_IND = N_IND,
                L_THIN = L_THIN, seeds = seeds)),
  "sbc_hierarchical_results.rds")
cat("Results saved.\n=== DONE ===\n")
