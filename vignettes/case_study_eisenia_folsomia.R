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

