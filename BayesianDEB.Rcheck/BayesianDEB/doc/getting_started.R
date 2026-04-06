## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----install------------------------------------------------------------------
# # Install cmdstanr (required backend)
# install.packages("cmdstanr",
#   repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()
# 
# # Install BayesianDEB
# # remotes::install_github("sciom/BayesianDEB")

## ----data---------------------------------------------------------------------
# library(BayesianDEB)
# data(eisenia_growth)
# 
# # Use first individual
# df1 <- eisenia_growth[eisenia_growth$id == 1, ]
# dat <- bdeb_data(growth = df1, f_food = 1.0)
# dat

## ----model--------------------------------------------------------------------
# mod <- bdeb_model(
#   data   = dat,
#   type   = "individual",
#   priors = list(
#     p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
#     p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
#     kappa   = prior_beta(a = 3, b = 2),
#     sigma_L = prior_halfnormal(sigma = 0.05)
#   )
# )
# mod

## ----fit----------------------------------------------------------------------
# fit <- bdeb_fit(mod, chains = 4, iter_sampling = 1000, seed = 123)
# fit

## ----diagnose-----------------------------------------------------------------
# bdeb_diagnose(fit)
# plot(fit, type = "trace")
# plot(fit, type = "pairs", pars = c("p_Am", "p_M", "kappa"))

## ----ppc----------------------------------------------------------------------
# ppc <- bdeb_ppc(fit, type = "growth")
# plot(ppc)

## ----derived------------------------------------------------------------------
# bdeb_derived(fit, quantities = c("L_inf", "growth_rate"))

## ----trajectory---------------------------------------------------------------
# plot(fit, type = "trajectory", n_draws = 200)

## ----hierarchical-------------------------------------------------------------
# dat_all <- bdeb_data(growth = eisenia_growth, f_food = 1.0)
# 
# mod_hier <- bdeb_model(dat_all, type = "hierarchical")
# 
# fit_hier <- bdeb_fit(mod_hier, chains = 4, adapt_delta = 0.9,
#                      threads_per_chain = 4)  # within-chain parallelism
# 
# bdeb_diagnose(fit_hier)
# bdeb_summary(fit_hier, pars = c("mu_log_p_Am", "sigma_log_p_Am", "p_M", "kappa"))

## ----debtox-------------------------------------------------------------------
# data(debtox_growth)
# 
# # Concentration mapping
# conc_map <- setNames(
#   c(0, 20, 80, 200),
#   c("1", "2", "3", "4")
# )
# 
# dat_tox <- bdeb_data(
#   growth = debtox_growth,
#   concentration = conc_map,
#   f_food = 1.0
# )
# 
# mod_tox <- bdeb_tox(dat_tox, stress = "assimilation")
# fit_tox <- bdeb_fit(mod_tox, chains = 4, adapt_delta = 0.95)
# 
# # EC50 and NEC
# bdeb_ec50(fit_tox)
# 
# # Dose-response plot
# plot_dose_response(fit_tox)

## ----priors-------------------------------------------------------------------
# # View defaults
# prior_default("individual")
# 
# # Override specific priors
# my_priors <- list(
#   p_Am  = prior_lognormal(mu = 2.0, sigma = 0.3),
#   kappa = prior_beta(a = 5, b = 2)
# )

## ----obs----------------------------------------------------------------------
# # Robust to outliers
# mod <- bdeb_model(dat, type = "individual",
#   observation = list(growth = obs_student_t(nu = 5)))
# 
# # Multiplicative error
# mod <- bdeb_model(dat, type = "individual",
#   observation = list(growth = obs_lognormal()))

