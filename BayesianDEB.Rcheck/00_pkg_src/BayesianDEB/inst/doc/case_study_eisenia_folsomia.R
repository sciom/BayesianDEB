## ----setup, include = FALSE---------------------------------------------------
HAS_CMDSTAN <- requireNamespace("cmdstanr", quietly = TRUE) &&
  isTRUE(nzchar(tryCatch(cmdstanr::cmdstan_path(),
                         error = function(e) "")))
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  eval = HAS_CMDSTAN
)
library(BayesianDEB)

## ----illustrative-note, eval = HAS_CMDSTAN, echo = FALSE, results = "asis"----
cat("> **Illustrative settings.** To keep build time within CRAN's",
    "vignette limits the chunks below use 2 chains × (300 + 300)",
    "iterations on a 5-individual subset of `eisenia_growth`.",
    "The publication-grade analysis with 4 chains × (1000 + 1000)",
    "iterations on all 21 individuals is reproduced by the scripts in the",
    "replication archive (`BayesianDEB_replication.zip`).\n\n")

## ----asis-warning, eval = !HAS_CMDSTAN, echo = FALSE, results = "asis"--------
# cat("> **Note.** This vignette requires `cmdstanr` and a working CmdStan",
#     "installation; the chunks below are not evaluated in the current",
#     "environment.  See `vignette(\"getting_started\", \"BayesianDEB\")`",
#     "for installation instructions.\n")

## ----prerequisites------------------------------------------------------------
library(BayesianDEB)
library(ggplot2)
library(posterior)  # for summarise_draws()

## ----check-stan---------------------------------------------------------------
# Internal helper; emits an informative error when CmdStan is missing.
BayesianDEB:::check_cmdstanr()

## ----eisenia-explore----------------------------------------------------------
data(eisenia_growth)

# Structure: 273 obs, 3 variables (id, time, length)
str(eisenia_growth)

length(unique(eisenia_growth$id))   # 21 individuals
length(unique(eisenia_growth$time)) # 13 time points (days 0–84)

## ----eisenia-plot, fig.cap="Growth trajectories of 21 *E. fetida* individuals.  Structural length $L = V^{1/3}$ measured weekly over 12 weeks."----
ggplot(eisenia_growth, aes(time, length, group = id)) +
  geom_line(alpha = 0.3, colour = "steelblue") +
  geom_point(size = 0.8, alpha = 0.4) +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = expression(paste("Structural length ", L, " (cm)")),
       title = "Eisenia fetida: 21 individuals, 12 weeks")

## ----ind-data-----------------------------------------------------------------
df1 <- eisenia_growth[eisenia_growth$id == 5, ]
dat1 <- bdeb_data(growth = df1, f_food = 1.0)
dat1

## ----ind-model----------------------------------------------------------------
mod1 <- bdeb_model(dat1, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    v       = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)
  ))
mod1

## ----ind-fit------------------------------------------------------------------
fit1 <- bdeb_fit(mod1,
  chains        = 2,
  iter_warmup   = 300,
  iter_sampling = 300,
  refresh       = 100,
  seed          = 42
)
fit1

## ----ind-diag-----------------------------------------------------------------
diag1 <- bdeb_diagnose(fit1)

## ----ind-trace, fig.cap="MCMC trace plots for core DEB parameters.  Well-mixed chains should appear as overlapping 'hairy caterpillars'."----
plot(fit1, type = "trace",
     pars = c("p_Am", "p_M", "kappa", "sigma_L"))

## ----ind-pairs, eval = HAS_CMDSTAN && requireNamespace("gridExtra", quietly = TRUE), fig.cap="Bivariate posterior scatter.  Strong correlation between $\\{p_{Am}\\}$ and $[p_M]$ is expected: both control ultimate size $L_\\infty = \\kappa \\{p_{Am}\\} / [p_M]$."----
# `bayesplot::mcmc_pairs` requires gridExtra (a Suggests of bayesplot).
plot(fit1, type = "pairs",
     pars = c("p_Am", "p_M", "kappa", "E_G"))

## ----ind-summary--------------------------------------------------------------
summary(fit1,
  pars = c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L"),
  prob = 0.95)

## ----ind-ppc, fig.cap="Posterior predictive check: grey lines are replicated growth trajectories, red points are observed data."----
ppc1 <- bdeb_ppc(fit1, type = "growth")
plot(ppc1, n_draws = 200)

## ----ind-traj, fig.cap="Posterior predicted trajectories (blue) with observed data (black points).  The spread reflects parameter uncertainty."----
plot(fit1, type = "trajectory", n_draws = 200)

## ----ind-derived--------------------------------------------------------------
der1 <- bdeb_derived(fit1,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"), f = 1.0)

summarise_draws(der1,
  "mean", "sd",
  "q2.5"  = ~quantile(.x, 0.025),
  "q97.5" = ~quantile(.x, 0.975))

## ----ind-food, fig.cap="Posterior distributions of $L_\\infty$ at $f = 1.0$ (blue) and $f = 0.7$ (orange)."----
d_f10 <- bdeb_derived(fit1, quantities = "L_inf", f = 1.0)
d_f07 <- bdeb_derived(fit1, quantities = "L_inf", f = 0.7)

df_compare <- data.frame(
  L_inf = c(d_f10$L_inf, d_f07$L_inf),
  food  = rep(c("f = 1.0", "f = 0.7"), each = nrow(d_f10))
)

ggplot(df_compare, aes(x = L_inf, fill = food)) +
  geom_density(alpha = 0.4) +
  theme_bw(base_size = 12) +
  labs(x = expression(L[infinity] ~ "(cm)"),
       y = "Posterior density",
       fill = "Food level")

## ----hier-data----------------------------------------------------------------
# Illustrative subset of 5 individuals; replication archive uses all 21.
dat_all <- bdeb_data(
  growth = eisenia_growth[eisenia_growth$id %in% 1:5, ],
  f_food = 1.0
)
dat_all

## ----hier-model---------------------------------------------------------------
mod_h <- bdeb_model(dat_all, type = "hierarchical",
  priors = list(
    mu_log_p_Am    = prior_normal(mu = 1.5, sigma = 0.5),
    sigma_log_p_Am = prior_exponential(rate = 2),
    p_M            = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa          = prior_beta(a = 3, b = 2),
    v              = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G            = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L        = prior_halfnormal(sigma = 0.05)
  ))

## ----hier-fit-----------------------------------------------------------------
fit_h <- bdeb_fit(mod_h,
  chains        = 2,
  iter_warmup   = 300,
  iter_sampling = 300,
  refresh       = 100,
  seed          = 123
)

## ----hier-diag----------------------------------------------------------------
bdeb_diagnose(fit_h)

## ----hier-trace, fig.cap="Trace plots for population-level hyperparameters $\\mu_{\\log p_{Am}}$ and $\\sigma_{\\log p_{Am}}$."----
plot(fit_h, type = "trace",
     pars = c("mu_log_p_Am", "sigma_log_p_Am"))

## ----hier-post, fig.cap="Marginal posterior densities for shared parameters."----
plot(fit_h, type = "posterior",
     pars = c("mu_log_p_Am", "sigma_log_p_Am", "p_M", "kappa"))

## ----hier-pop-----------------------------------------------------------------
summary(fit_h,
  pars = c("mu_log_p_Am", "sigma_log_p_Am",
           "p_M", "kappa", "v", "E_G", "sigma_L"),
  prob = 0.95)

## ----hier-shrinkage, fig.cap="Individual-level $\\{p_{Am}\\}$ estimates (points: posterior means; bars: 90% CI) compared to the population mean (dashed red line).  Shrinkage toward the mean is visible for individuals with noisier data."----
n_ind <- dat_all$n_ind
ind_summary <- summary(fit_h,
  pars = paste0("p_Am_ind[", seq_len(n_ind), "]"),
  prob = 0.90)

pop_summary <- summary(fit_h, pars = "mu_log_p_Am")
pop_mean_pAm <- exp(as.data.frame(pop_summary)$mean)

ind_df <- as.data.frame(ind_summary)
ind_df$individual <- seq_len(n_ind)

ggplot(ind_df, aes(x = individual, y = mean)) +
  geom_pointrange(aes(ymin = `5%`, ymax = `95%`),
                  colour = "steelblue", size = 0.4) +
  geom_hline(yintercept = pop_mean_pAm, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(x = "Individual", y = expression({p[Am]} ~ "(J/d/cm"^2*")"),
       title = "Individual assimilation rates with 90% CI")

## ----hier-new-----------------------------------------------------------------
summary(fit_h, pars = "p_Am_new", prob = 0.95)

## ----compare------------------------------------------------------------------
s_ind  <- summary(fit1, pars = c("p_Am", "p_M", "kappa"), prob = 0.90)
s_hier <- summary(fit_h, pars = c("p_M", "kappa"), prob = 0.90)

cat("=== Individual model (id = 5, n = 1) ===\n")
print(as.data.frame(s_ind), digits = 3, row.names = FALSE)

cat("\n=== Hierarchical model (n = 21) ===\n")
print(as.data.frame(s_hier), digits = 3, row.names = FALSE)

## ----debtox-explore, fig.cap="Growth trajectories under 4 toxicant concentrations.  Higher concentrations suppress growth through reduced assimilation."----
data(debtox_growth)

ggplot(debtox_growth,
       aes(time, length, colour = factor(concentration), group = id)) +
  geom_line(alpha = 0.3) +
  geom_point(size = 0.8, alpha = 0.4) +
  facet_wrap(~concentration, labeller = label_both) +
  theme_bw(base_size = 11) +
  scale_colour_brewer(palette = "RdYlBu", direction = -1) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       colour = "Concentration") +
  theme(legend.position = "none")

## ----debtox-prep--------------------------------------------------------------
conc_levels <- unique(debtox_growth$concentration)
conc_map <- setNames(conc_levels, as.character(conc_levels))

dat_tox <- bdeb_data(
  growth        = debtox_growth,
  concentration = conc_map,
  f_food        = 1.0
)
dat_tox

## ----debtox-model-------------------------------------------------------------
mod_tox <- bdeb_tox(dat_tox, stress = "assimilation",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    v       = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05),
    k_d     = prior_lognormal(mu = -1.0, sigma = 1.0),
    z_w     = prior_lognormal(mu = 2.5, sigma = 1.0),
    b_w     = prior_lognormal(mu = -5.0, sigma = 2.0)
  ))
mod_tox

## ----debtox-fit---------------------------------------------------------------
fit_tox <- bdeb_fit(mod_tox, algorithm = "variational",
                    seed = 77, refresh = 0)

## ----debtox-trace-tox, eval = FALSE, fig.cap="Trace plots for the three toxicological parameters.  Good mixing is essential for reliable EC$_{50}$ and NEC estimates."----
# # Re-fit with algorithm = "sampling" for these diagnostics; see
# # the replication archive for the publication-grade analysis.
# plot(fit_tox, type = "trace", pars = c("k_d", "z_w", "b_w"))

## ----debtox-post-tox, fig.cap="Marginal posterior densities for toxicological parameters."----
plot(fit_tox, type = "posterior", pars = c("k_d", "z_w", "b_w"))

## ----debtox-pairs, eval = HAS_CMDSTAN && requireNamespace("gridExtra", quietly = TRUE), fig.cap="Posterior pairs for toxicological parameters.  A correlation between $z_w$ and $b_w$ is expected since both determine the shape of the dose-response curve."----
# `bayesplot::mcmc_pairs` requires gridExtra (a Suggests of bayesplot).
plot(fit_tox, type = "pairs", pars = c("z_w", "b_w", "k_d"))

## ----debtox-summary-----------------------------------------------------------
summary(fit_tox,
  pars = c("p_Am", "p_M", "kappa", "v", "E_G",
           "k_d", "z_w", "b_w", "sigma_L"),
  prob = 0.95)

## ----debtox-ec50, fig.cap="Posterior distribution of EC$_{50}$ (blue histogram) with the posterior median (red dashed line).  The full distribution — not just a point estimate — is available for regulatory risk assessment."----
ec <- bdeb_ec50(fit_tox, prob = 0.95)
print(ec$summary, digits = 3)

hist(ec$draws, breaks = 50, col = "steelblue", border = "white",
     main = expression("Posterior distribution of EC"[50]),
     xlab = "Concentration", freq = FALSE)
abline(v = ec$summary$median[1], col = "red", lwd = 2, lty = 2)
legend("topright", "Posterior median",
       col = "red", lty = 2, lwd = 2, bty = "n")

## ----debtox-dr, fig.cap="Dose-response curve with posterior uncertainty bands (blue lines: individual posterior draws).  The dashed horizontal line marks 50% effect; vertical dashed lines mark the NEC (green) and EC$_{50}$ (red).  The vignette uses lite settings; replication archive uses 200 draws."----
plot_dose_response(fit_tox, n_draws = 30, n_conc = 25)

## ----debtox-sens--------------------------------------------------------------
mod_tox2 <- bdeb_tox(dat_tox, stress = "assimilation",
  priors = list(
    z_w = prior_lognormal(mu = 3.0, sigma = 0.3),  # tighter
    b_w = prior_lognormal(mu = -5.0, sigma = 2.0)
  ))
fit_tox2 <- bdeb_fit(mod_tox2, algorithm = "variational",
                     seed = 78, refresh = 0)

cat("=== Original: z_w ~ LogNormal(2.5, 1.0) ===\n")
summary(fit_tox,  pars = c("z_w", "b_w"), prob = 0.95)

cat("\n=== Tighter:  z_w ~ LogNormal(3.0, 0.3) ===\n")
summary(fit_tox2, pars = c("z_w", "b_w"), prob = 0.95)

## ----convert-length-----------------------------------------------------------
# Example: measured body lengths in mm for E. fetida
L_physical_mm <- c(12, 18, 25, 30)
delta_M <- 0.24

# Convert to structural length in cm
L_structural_cm <- delta_M * L_physical_mm / 10
L_structural_cm
# [1] 0.288 0.432 0.600 0.720

## ----prior-pred---------------------------------------------------------------
set.seed(42)
n_sim <- 4000

# Sample from priors
p_Am_sim  <- rlnorm(n_sim, 1.5, 0.5)
p_M_sim   <- rlnorm(n_sim, -1.0, 0.5)
kappa_sim <- rbeta(n_sim, 3, 2)
v_sim     <- rlnorm(n_sim, -1.5, 0.5)
E_G_sim   <- rlnorm(n_sim, 6.0, 0.5)

# Prior predictive for L_inf
L_inf_prior <- kappa_sim * p_Am_sim / p_M_sim

hist(L_inf_prior, breaks = 50, col = "steelblue", border = "white",
     main = "Prior predictive: ultimate structural length",
     xlab = expression(L[infinity] ~ "(cm)"), xlim = c(0, 50))
# Should cover plausible range for earthworms (~2-20 cm structural)

## ----obs-models---------------------------------------------------------------
# Robust to outliers: Student-t with 5 df
mod_robust <- bdeb_model(dat1, type = "individual",
  observation = list(growth = obs_student_t(nu = 5)))

# Multiplicative error (constant CV)
mod_logn <- bdeb_model(dat1, type = "individual",
  observation = list(growth = obs_lognormal()))

## ----arrhenius----------------------------------------------------------------
# Experiment at 22 C, reference 20 C, typical T_A for ectotherms
cT <- arrhenius(temp = 273.15 + 22, T_ref = 273.15 + 20, T_A = 8000)
cat("Temperature correction factor:", round(cT, 3), "\n")
# Rate at reference temperature: p_Am_ref = p_Am_obs / cT

## ----fluxes-------------------------------------------------------------------
fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                 p_Am = 5, p_M = 0.5, kappa = 0.75,
                 v = 0.2, E_G = 400)

cat(sprintf("Assimilation  (p_A): %.3f J/d\n", fl$p_A))
cat(sprintf("Mobilisation  (p_C): %.3f J/d\n", fl$p_C))
cat(sprintf("Maintenance   (p_M): %.3f J/d\n", fl$p_M))
cat(sprintf("Growth        (p_G): %.3f J/d\n", fl$p_G))
cat(sprintf("Struct. length (L) : %.3f cm\n",  fl$L))
cat(sprintf("Scaled reserve (e) : %.3f\n",     fl$e))

## ----repro-convert------------------------------------------------------------
cumul <- data.frame(
  id = rep(1, 5),
  time = c(0, 7, 14, 21, 28),
  cumulative = c(0, 10, 30, 60, 100)
)
repro_to_intervals(cumul)
#   id t_start t_end count
# 1  1       0     7    10
# 2  1       7    14    20
# 3  1      14    21    30
# 4  1      21    28    40

## ----session------------------------------------------------------------------
sessionInfo()

