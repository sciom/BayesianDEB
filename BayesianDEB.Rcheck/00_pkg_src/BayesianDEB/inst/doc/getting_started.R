## ----setup, include = FALSE---------------------------------------------------
HAS_CMDSTAN <- requireNamespace("cmdstanr", quietly = TRUE) &&
  isTRUE(nzchar(tryCatch(cmdstanr::cmdstan_path(),
                         error = function(e) "")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4,
  eval = HAS_CMDSTAN
)
library(BayesianDEB)

## ----install, eval = FALSE----------------------------------------------------
# # Install cmdstanr (required backend)
# install.packages("cmdstanr",
#   repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()
# 
# # Install BayesianDEB
# # remotes::install_github("sciom/BayesianDEB")

## ----data---------------------------------------------------------------------
data(eisenia_growth)

# Use first individual
df1 <- eisenia_growth[eisenia_growth$id == 1, ]
dat <- bdeb_data(growth = df1, f_food = 1.0)
dat

## ----model--------------------------------------------------------------------
mod <- bdeb_model(
  data   = dat,
  type   = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)
  )
)
mod

## ----fit----------------------------------------------------------------------
fit <- bdeb_fit(mod,
                chains = 2, iter_warmup = 300, iter_sampling = 300,
                seed = 123, refresh = 100)
fit

## ----diagnose-----------------------------------------------------------------
bdeb_diagnose(fit)
plot(fit, type = "trace")

## ----diagnose-pairs, eval = HAS_CMDSTAN && requireNamespace("gridExtra", quietly = TRUE)----
# `bayesplot::mcmc_pairs` requires gridExtra (a Suggests of bayesplot).
plot(fit, type = "pairs", pars = c("p_Am", "p_M", "kappa"))

## ----ppc----------------------------------------------------------------------
ppc <- bdeb_ppc(fit, type = "growth")
plot(ppc)

## ----derived------------------------------------------------------------------
bdeb_derived(fit, quantities = c("L_inf", "growth_rate"))

## ----trajectory---------------------------------------------------------------
plot(fit, type = "trajectory", n_draws = 200)

## ----hierarchical-------------------------------------------------------------
dat_all <- bdeb_data(
  growth = eisenia_growth[eisenia_growth$id %in% 1:5, ],
  f_food = 1.0
)

mod_hier <- bdeb_model(dat_all, type = "hierarchical")

fit_hier <- bdeb_fit(mod_hier,
                     chains = 2, iter_warmup = 300, iter_sampling = 300,
                     seed = 123, refresh = 100)

bdeb_diagnose(fit_hier)
summary(fit_hier, pars = c("mu_log_p_Am", "sigma_log_p_Am", "p_M", "kappa"))

## ----debtox-data--------------------------------------------------------------
data(debtox_growth)

# Concentration mapping
conc_map <- setNames(
  c(0, 20, 80, 200),
  c("1", "2", "3", "4")
)

dat_tox <- bdeb_data(
  growth = debtox_growth,
  concentration = conc_map,
  f_food = 1.0
)

mod_tox <- bdeb_tox(dat_tox, stress = "assimilation")

## ----debtox-fit---------------------------------------------------------------
fit_tox <- bdeb_fit(mod_tox, algorithm = "variational",
                    seed = 123, refresh = 0)

## ----debtox-ec50--------------------------------------------------------------
bdeb_ec50(fit_tox)

## ----debtox-plot--------------------------------------------------------------
plot_dose_response(fit_tox, n_draws = 20, n_conc = 25, dt = 1.0)

## ----priors-------------------------------------------------------------------
# View defaults
prior_default("individual")

# Override specific priors
my_priors <- list(
  p_Am  = prior_lognormal(mu = 2.0, sigma = 0.3),
  kappa = prior_beta(a = 5, b = 2)
)

## ----obs----------------------------------------------------------------------
# Robust to outliers
mod <- bdeb_model(dat, type = "individual",
  observation = list(growth = obs_student_t(nu = 5)))

# Multiplicative error
mod <- bdeb_model(dat, type = "individual",
  observation = list(growth = obs_lognormal()))

