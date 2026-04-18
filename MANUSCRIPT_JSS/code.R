## =============================================================
## BayesianDEB: Replication script for JSS manuscript
## Hackenberger, Djerdj & Hackenberger (2026)
##
## This script reproduces all tables and figures in the paper.
## Requires: BayesianDEB (>= 0.1.2), cmdstanr, CmdStan,
##           posterior, ggplot2, GGally, deSolve, tidyr
## Runtime: ~3 hours on a 16-core machine
##
## Long-running SBC scripts are separate (see below).
## =============================================================

options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)

library("BayesianDEB")
library("posterior")
library("ggplot2")

set.seed(42)

## --- Consistent colour palette across all figures --------------------
col_posterior <- "#2166AC"   # blue for posterior draws / trajectories
col_observed  <- "#B2182B"   # red for observed data points
col_replicated <- "grey60"   # grey for PPC replicated data
col_prior     <- "#D6604D"   # light red for prior densities


## =============================================================
## Section 5.1: Individual growth — simulated E. fetida
## =============================================================

cat("=== Section 5.1: Individual growth (simulated) ===\n")

data("eisenia_growth")
df1 <- eisenia_growth[eisenia_growth$id == 5, ]
dat <- bdeb_data(growth = df1, f_food = 1.0)

mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)))

fit <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
  adapt_delta = 0.9, seed = 42, refresh = 200)

bdeb_diagnose(fit)

cat("\n--- Derived quantities ---\n")
der <- bdeb_derived(fit,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"),
  f = 1.0)
print(as.data.frame(summarise_draws(der, "mean", "sd",
  "q5" = ~quantile(.x, 0.05),
  "q95" = ~quantile(.x, 0.95))), digits = 3)

## LOO comparison (Gaussian vs. log-normal)
mod_ln <- bdeb_model(dat, type = "individual",
  observation = list(growth = obs_lognormal()),
  priors = list(
    p_Am  = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)))
fit_ln <- bdeb_fit(mod_ln, chains = 4,
  iter_sampling = 2000, adapt_delta = 0.9, seed = 42,
  refresh = 200)

loo_gauss <- bdeb_loo(fit)
loo_ln    <- bdeb_loo(fit_ln)
cat("\n--- LOO comparison ---\n")
print(loo::loo_compare(list(gaussian = loo_gauss,
  lognormal = loo_ln)))


## =============================================================
## Section 5.2: Growth with reproduction — simulated F. candida
## =============================================================

cat("\n=== Section 5.2: Growth with reproduction ===\n")

growth_fc <- data.frame(id = 1,
  time = seq(0, 28, by = 2),
  length = c(0.048, 0.053, 0.058, 0.063, 0.067, 0.071,
             0.075, 0.078, 0.081, 0.084, 0.086, 0.088,
             0.089, 0.090, 0.091))
repro_fc <- data.frame(id = 1, t_start = 0,
  t_end = 28, count = 118)

dat_gr <- bdeb_data(growth = growth_fc,
  reproduction = repro_fc, f_food = 1.0)

mod_gr <- bdeb_model(dat_gr, type = "growth_repro",
  priors = prior_species("Folsomia_candida",
    type = "growth_repro"),
  observation = list(growth = obs_normal(),
    reproduction = obs_negbinom()))

fit_gr <- bdeb_fit(mod_gr, chains = 4,
  iter_sampling = 2000, adapt_delta = 0.9, seed = 42,
  refresh = 200)

cat("\n--- Growth+repro diagnostics ---\n")
bdeb_diagnose(fit_gr)

cat("\n--- Growth+repro posterior ---\n")
print(as.data.frame(bdeb_summary(fit_gr,
  pars = c("p_Am", "p_M", "kappa", "k_J", "k_R", "phi_R",
           "sigma_L"), prob = 0.90)), digits = 3)


## =============================================================
## Section 5.3: Hierarchical model
## =============================================================

cat("\n=== Section 5.3: Hierarchical model ===\n")

dat_all <- bdeb_data(growth = eisenia_growth, f_food = 1.0)

mod_h <- bdeb_model(dat_all, type = "hierarchical",
  priors = list(
    mu_log_p_Am    = prior_normal(mu = 1.5, sigma = 0.5),
    sigma_log_p_Am = prior_exponential(rate = 2)))

fit_h <- bdeb_fit(mod_h, chains = 4, adapt_delta = 0.95,
  threads_per_chain = 4, seed = 123, refresh = 200)

cat("\n--- Hierarchical diagnostics ---\n")
bdeb_diagnose(fit_h)


## =============================================================
## Section 5.4: DEBtox — Van Gestel (1991) Cd data
## =============================================================

cat("\n=== Section 5.4: DEBtox — Van Gestel 1991 ===\n")

data("eisenia_cd")
conc_map <- setNames(unique(eisenia_cd$concentration),
  as.character(unique(eisenia_cd$concentration)))
dat_cd <- bdeb_data(growth = eisenia_cd,
  concentration = conc_map, f_food = 1.0)

## Raw data figure
p_debtox_raw <- ggplot(eisenia_cd,
  aes(time, length, colour = factor(concentration))) +
  geom_point(size = 2) + geom_line() +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       colour = "Cd (mg/kg)")
ggsave("fig_debtox_rawdata.pdf", p_debtox_raw,
  width = 7, height = 4)

mod_cd <- bdeb_tox(dat_cd, stress = "assimilation",
  priors = list(
    z_w = prior_lognormal(mu = 3.0, sigma = 1.0),
    b_w = prior_lognormal(mu = -4.0, sigma = 1.5)))

fit_cd <- bdeb_fit(mod_cd, chains = 4, iter_warmup = 1000,
  iter_sampling = 2000, adapt_delta = 0.95,
  max_treedepth = 12, seed = 42, refresh = 200)

## Table 4: DEBtox posterior (tab:debtox)
cat("\n--- Table 4: DEBtox posterior ---\n")
print(as.data.frame(bdeb_summary(fit_cd,
  pars = c("p_Am","p_M","kappa","k_d","z_w","b_w",
           "sigma_L"), prob = 0.90)), digits = 3)

ec <- bdeb_ec50(fit_cd, prob = 0.90)
cat("\n--- EC50/NEC ---\n")
print(ec$summary, digits = 3)

## Dose-response figure
p_dr <- plot_dose_response(fit_cd, n_draws = 100, seed = 1)
ggsave("fig_debtox_doseresponse.pdf", p_dr,
  width = 6, height = 4)


## =============================================================
## Section 6.1: SBC
## =============================================================

cat("\n=== Section 6.1: SBC ===\n")
cat("SBC scripts are run separately (long-running):\n")
cat("  sbc_validation.R     — individual model, 500 reps, ~90 min\n")
cat("  sbc_hierarchical.R   — hierarchical model, 50 reps, ~5 hours\n")

if (file.exists("sbc_results.rds")) {
  sbc <- readRDS("sbc_results.rds")
  cat("SBC individual: loaded from sbc_results.rds\n")
  cat(sprintf("  Valid: %d / %d\n", sbc$n_valid, sbc$n_total))
}
if (file.exists("sbc_hierarchical_results.rds")) {
  sbc_h <- readRDS("sbc_hierarchical_results.rds")
  cat("SBC hierarchical: loaded from sbc_hierarchical_results.rds\n")
}


## =============================================================
## Section 6.2: Identifiability / contraction (Table 5)
## =============================================================

cat("\n=== Section 6.2: Identifiability analysis ===\n")
cat("Running identifiability_analysis.R...\n")
source("identifiability_analysis.R")


## =============================================================
## Section 6.3: Real data — Neuhauser (1980)
## =============================================================

cat("\n=== Section 6.3: Neuhauser 1980 real data ===\n")

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
dat_real <- bdeb_data(growth = eisenia_real, f_food = 1.0)
mod_real <- bdeb_model(dat_real, type = "individual",
  priors = list(
    p_Am  = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M   = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa = prior_beta(a = 3, b = 2),
    v     = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G   = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)),
  temperature = list(T = 298.15, T_ref = 293.15, T_A = 8000))

fit_real <- bdeb_fit(mod_real, chains = 4, iter_sampling = 2000,
  adapt_delta = 0.95, max_treedepth = 12,
  seed = 42, refresh = 200)

## Table 6: Real data posterior (tab:realdata)
cat("\n--- Table 6: Real data posterior ---\n")
s_real <- bdeb_summary(fit_real,
  pars = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"),
  prob = 0.90)
print(as.data.frame(s_real), digits = 3)

## Derived quantities
cat("\n--- Derived quantities ---\n")
der_real <- bdeb_derived(fit_real,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"),
  f = 1.0)
print(as.data.frame(summarise_draws(der_real, "mean", "sd",
  "q5" = ~quantile(.x, 0.05),
  "q95" = ~quantile(.x, 0.95))), digits = 3)

## Table: posterior correlation matrix (tab:corr)
cat("\n--- Posterior correlation (Spearman) ---\n")
draws_mat <- as_draws_matrix(fit_real$fit$draws())
pars_draws <- draws_mat[, c("p_Am","p_M","kappa","v","E_G","sigma_L")]
print(round(cor(pars_draws, method = "spearman"), 2))

## Figures: diagnostics
draws_df <- as_draws_df(fit_real$fit$draws())
draws_clean <- draws_df[complete.cases(
  draws_df[, c("p_Am","p_M","kappa","v","E_G","sigma_L")]), ]

## Trace plots (fig:diagnostics)
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

## Posterior densities (fig:posterior)
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

## Pairs plot (fig:pairs)
pairs_df <- as.data.frame(
  draws_clean[, c("p_Am","p_M","kappa","E_G")])
p_pairs <- GGally::ggpairs(pairs_df,
  upper = list(continuous = "cor"),
  lower = list(continuous = "points"),
  diag = list(continuous = "densityDiag")) +
  theme_bw(base_size = 10)
ggsave("fig_pairs.pdf", p_pairs, width = 8, height = 8)

## Trajectory (fig:trajectory_real)
L_hat_vars <- grep("^L_hat", names(draws_clean), value = TRUE)
L_hat <- as.matrix(draws_clean[, L_hat_vars])
L_hat <- L_hat[complete.cases(L_hat), ]
idx <- sort(sample.int(nrow(L_hat), min(200, nrow(L_hat))))
t_obs <- dat_real$growth$time
traj_list <- lapply(idx, function(i)
  data.frame(time = t_obs[1:ncol(L_hat)],
             length = L_hat[i,], draw = i))
traj_df <- do.call(rbind, traj_list)
obs_df <- data.frame(time = t_obs, length = dat_real$growth$length)

p_traj <- ggplot() +
  geom_line(data = traj_df,
    aes(time, length, group = draw),
    alpha = 0.1, colour = col_posterior) +
  geom_point(data = obs_df, aes(time, length), size = 2.5) +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "Posterior Predicted Trajectories")
ggsave("fig_trajectory.pdf", p_traj, width = 6, height = 4)

## PPC (fig:ppc)
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
    alpha = 0.1, colour = col_replicated) +
  geom_point(data = obs_df, aes(time, length),
    size = 2.5, colour = col_observed) +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "Posterior Predictive Check")
ggsave("fig_ppc.pdf", p_ppc, width = 6, height = 4)

## Prior predictive (fig:prior_pred)
pp <- bdeb_prior_predictive(mod_real, n_draws = 500, seed = 42)
p_prior <- plot(pp, n_draws = 100)
ggsave("fig_prior_predictive.pdf", p_prior, width = 6, height = 4)

## Prior vs posterior (fig:prior_post)
p_pp <- plot(fit_real, type = "prior_posterior")
ggsave("fig_prior_posterior.pdf", p_pp, width = 8, height = 5)


## =============================================================
## Section 6.4: Prior sensitivity analysis (Table 7)
## =============================================================

cat("\n=== Section 6.4: Prior sensitivity analysis ===\n")

## Prior A: informative (sigma = 0.5)
dat_sens <- bdeb_data(
  growth = eisenia_growth[eisenia_growth$id == 5, ],
  f_food = 1.0)

mod_A <- bdeb_model(dat_sens, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    v       = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)))

fit_A <- bdeb_fit(mod_A, chains = 4, iter_sampling = 2000,
  adapt_delta = 0.9, seed = 42, refresh = 200)

## Prior B: weakly informative (sigma = 1.0, defaults)
mod_B <- bdeb_model(dat_sens, type = "individual",
  priors = list(
    sigma_L = prior_halfnormal(sigma = 0.1)))

fit_B <- bdeb_fit(mod_B, chains = 4, iter_sampling = 2000,
  adapt_delta = 0.9, seed = 42, refresh = 200)

cat("\n--- Table 7: Prior sensitivity ---\n")
cat("Prior A (informative, sigma = 0.5):\n")
print(as.data.frame(bdeb_summary(fit_A,
  pars = c("p_Am","p_M","kappa","v","sigma_L"),
  prob = 0.90)), digits = 3)

cat("\nPrior B (weakly informative, sigma = 1.0):\n")
print(as.data.frame(bdeb_summary(fit_B,
  pars = c("p_Am","p_M","kappa","v","sigma_L"),
  prob = 0.90)), digits = 3)

## Derived quantities for both
draws_A <- as.data.frame(as_draws_matrix(fit_A$fit$draws()))
draws_B <- as.data.frame(as_draws_matrix(fit_B$fit$draws()))

cat("\nDerived L_inf — Prior A:\n")
L_inf_A <- draws_A$kappa * draws_A$p_Am / draws_A$p_M
cat(sprintf("  median = %.2f, 90%% CI = [%.2f, %.2f]\n",
  median(L_inf_A),
  quantile(L_inf_A, 0.05), quantile(L_inf_A, 0.95)))

cat("Derived L_inf — Prior B:\n")
L_inf_B <- draws_B$kappa * draws_B$p_Am / draws_B$p_M
cat(sprintf("  median = %.2f, 90%% CI = [%.2f, %.2f]\n",
  median(L_inf_B),
  quantile(L_inf_B, 0.05), quantile(L_inf_B, 0.95)))


## =============================================================
## Session info
## =============================================================

cat("\n=== Session info ===\n")
bdeb_session_info(fit_real)
cat("\n")
sessionInfo()
