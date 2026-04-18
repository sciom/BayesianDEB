## =============================================================
## BayesianDEB: Replication script for JSS manuscript
## Hackenberger, Djerdj & Hackenberger (2026)
##
## This script reproduces all tables and figures in the paper.
## Requires: BayesianDEB (>= 0.1.1), cmdstanr, CmdStan,
##           posterior, ggplot2, deSolve, GGally
## Runtime: ~2 hours on a 16-core machine
## =============================================================

options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)

library("BayesianDEB")
library("posterior")
library("ggplot2")

set.seed(42)

## =============================================================
## Section 6.3: Real data case study — Neuhauser (1980)
## =============================================================

cat("=== Section 6.3: Neuhauser 1980 real data ===\n")

## Load bundled dataset (originally from EGrowth, Mathieu 2018)
data("eisenia_neuhauser")
eisenia_real <- eisenia_neuhauser

## Figure: raw data
p_raw <- ggplot(eisenia_real, aes(time, length)) +
  geom_point(size = 2.5) + geom_line(alpha = 0.3) +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "E. fetida growth (Neuhauser et al. 1980)")
ggsave("fig_rawdata.pdf", p_raw, width = 6, height = 4)

## Fit individual model with Arrhenius correction (25C -> 20C ref)
dat <- bdeb_data(growth = eisenia_real, f_food = 1.0)
mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am  = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M   = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa = prior_beta(a = 3, b = 2),
    v     = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G   = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)),
  temperature = list(T = 298.15, T_ref = 293.15, T_A = 8000))

fit <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
                adapt_delta = 0.95, max_treedepth = 12,
                seed = 42, refresh = 200)

## Table: posterior summary (Table realdata in paper)
cat("\n--- Table: Real data posterior ---\n")
s <- bdeb_summary(fit,
  pars = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"),
  prob = 0.90)
print(as.data.frame(s), digits = 3)

## Table: derived quantities
cat("\n--- Derived quantities ---\n")
der <- bdeb_derived(fit,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"),
  f = 1.0)
print(as.data.frame(summarise_draws(der, "mean", "sd",
  "q5" = ~quantile(.x, 0.05),
  "q95" = ~quantile(.x, 0.95))), digits = 3)

## Table: posterior correlation matrix (Table corr)
cat("\n--- Posterior correlation (Spearman) ---\n")
draws_mat <- as_draws_matrix(fit$fit$draws())
pars_draws <- draws_mat[, c("p_Am","p_M","kappa","v","E_G","sigma_L")]
print(round(cor(pars_draws, method = "spearman"), 2))

## Figures: diagnostics
draws_df <- as_draws_df(fit$fit$draws())
draws_clean <- draws_df[complete.cases(
  draws_df[, c("p_Am","p_M","kappa","v","E_G","sigma_L")]), ]

# Trace
pars_trace <- c("p_Am", "p_M", "kappa", "sigma_L")
trace_long <- tidyr::pivot_longer(
  draws_clean[, c(pars_trace, ".chain", ".iteration")],
  cols = all_of(pars_trace),
  names_to = "parameter", values_to = "value")
p_trace <- ggplot(trace_long,
  aes(.iteration, value, colour = factor(.chain))) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~parameter, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(x = "Iteration", y = "Value", colour = "Chain") +
  theme(legend.position = "bottom")
ggsave("fig_trace.pdf", p_trace, width = 8, height = 5)

# Posterior densities
pars_all <- c("p_Am","p_M","kappa","v","E_G","sigma_L")
dens_long <- tidyr::pivot_longer(
  draws_clean[, c(pars_all, ".chain")],
  cols = all_of(pars_all),
  names_to = "parameter", values_to = "value")
p_post <- ggplot(dens_long, aes(value, fill = factor(.chain))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw(base_size = 11) +
  labs(x = "Value", y = "Density", fill = "Chain") +
  theme(legend.position = "bottom")
ggsave("fig_posterior.pdf", p_post, width = 8, height = 5)

# Pairs
pairs_df <- as.data.frame(
  draws_clean[, c("p_Am","p_M","kappa","E_G")])
p_pairs <- GGally::ggpairs(pairs_df,
  upper = list(continuous = "cor"),
  lower = list(continuous = "points"),
  diag = list(continuous = "densityDiag")) +
  theme_bw(base_size = 10)
ggsave("fig_pairs.pdf", p_pairs, width = 8, height = 8)

# Trajectory
L_hat_vars <- grep("^L_hat", names(draws_clean), value = TRUE)
L_hat <- as.matrix(draws_clean[, L_hat_vars])
L_hat <- L_hat[complete.cases(L_hat), ]
idx <- sort(sample.int(nrow(L_hat), min(200, nrow(L_hat))))
t_obs <- dat$growth$time
traj_list <- lapply(idx, function(i)
  data.frame(time = t_obs[1:ncol(L_hat)],
             length = L_hat[i,], draw = i))
traj_df <- do.call(rbind, traj_list)
obs_df <- data.frame(time = t_obs, length = dat$growth$length)

p_traj <- ggplot() +
  geom_line(data = traj_df,
    aes(time, length, group = draw),
    alpha = 0.1, colour = "steelblue") +
  geom_point(data = obs_df, aes(time, length), size = 2.5) +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "Posterior Predicted Trajectories")
ggsave("fig_trajectory.pdf", p_traj, width = 6, height = 4)

# PPC
L_rep_vars <- grep("^L_rep", names(draws_clean), value = TRUE)
L_rep <- as.matrix(draws_clean[, L_rep_vars])
L_rep <- L_rep[complete.cases(L_rep), ]
idx_r <- sort(sample.int(nrow(L_rep), min(200, nrow(L_rep))))
ppc_list <- lapply(idx_r, function(i)
  data.frame(time = t_obs[1:ncol(L_rep)],
             length = L_rep[i,], draw = i))
ppc_df <- do.call(rbind, ppc_list)

p_ppc <- ggplot() +
  geom_line(data = ppc_df,
    aes(time, length, group = draw),
    alpha = 0.1, colour = "grey50") +
  geom_point(data = obs_df, aes(time, length),
    size = 2.5, colour = "red") +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "Posterior Predictive Check")
ggsave("fig_ppc.pdf", p_ppc, width = 6, height = 4)

# Prior predictive
pp <- bdeb_prior_predictive(mod, n_draws = 500, seed = 42)
p_prior <- plot(pp, n_draws = 100)
ggsave("fig_prior_predictive.pdf", p_prior, width = 6, height = 4)


## =============================================================
## Section 5.3: DEBtox — Van Gestel (1991)
## =============================================================

cat("\n=== Section 5.3: DEBtox — Van Gestel 1991 ===\n")

data("eisenia_cd")
conc_map <- setNames(unique(eisenia_cd$concentration),
  as.character(unique(eisenia_cd$concentration)))
dat_cd <- bdeb_data(growth = eisenia_cd,
  concentration = conc_map, f_food = 1.0)

# Raw data figure
p_debtox_raw <- ggplot(eisenia_cd,
  aes(time, length, colour = factor(concentration))) +
  geom_point(size = 2) + geom_line() +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       colour = "Cd (mg/kg)")
ggsave("fig_debtox_rawdata.pdf", p_debtox_raw, width = 7, height = 4)

mod_cd <- bdeb_tox(dat_cd, stress = "assimilation",
  priors = list(
    z_w = prior_lognormal(mu = 3.0, sigma = 1.0),
    b_w = prior_lognormal(mu = -4.0, sigma = 1.5)))

fit_cd <- bdeb_fit(mod_cd, chains = 4, iter_warmup = 1000,
  iter_sampling = 2000, adapt_delta = 0.95,
  max_treedepth = 12, seed = 42, refresh = 200)

cat("\n--- DEBtox posterior ---\n")
print(as.data.frame(bdeb_summary(fit_cd,
  pars = c("p_Am","p_M","kappa","k_d","z_w","b_w",
           "sigma_L"), prob = 0.90)), digits = 3)

ec <- bdeb_ec50(fit_cd, prob = 0.90)
cat("\n--- EC50/NEC ---\n")
print(ec$summary, digits = 3)

# Dose-response figure
p_dr <- plot_dose_response(fit_cd, n_draws = 100, seed = 1)
ggsave("fig_debtox_doseresponse.pdf", p_dr, width = 6, height = 4)


## =============================================================
## Section 6.1: SBC (individual model)
## =============================================================

cat("\n=== Section 6.1: SBC ===\n")
cat("SBC scripts: sbc_validation.R (individual, ~90 min)\n")
cat("             sbc_hierarchical.R (hierarchical, ~5 hours)\n")
cat("Pre-computed results in sbc_results.rds\n")

if (file.exists("sbc_results.rds")) {
  sbc <- readRDS("sbc_results.rds")
  cat("SBC individual: loaded from sbc_results.rds\n")
}


## =============================================================
## Session info
## =============================================================

cat("\n=== Session info ===\n")
bdeb_session_info(fit)
cat("\n")
sessionInfo()
