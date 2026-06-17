## *******************************************************************
## Quantitative comparison: BayesianDEB vs. deBInfer
## Fits the same DEB growth model to the same dataset using both
## Dataset: eisenia_growth (individual 5, 13 weekly observations)
## Runtime: cca 10 minutes
## *******************************************************************

library(BayesianDEB)
library(deBInfer)
library(deSolve)
library(posterior)

set.seed(42)

## *******************************************************************
## 1. BayesianDEB 
## *******************************************************************

cat("*** BayesianDEB fit ***\n")

data("eisenia_growth")
df1 <- eisenia_growth[eisenia_growth$id == 5, ]

dat <- bdeb_data(growth = df1, f_food = 1.0)
mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    v       = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)))

t_bdeb <- system.time({
  fit_bdeb <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
    adapt_delta = 0.9, seed = 42, refresh = 0)
})

draws_bdeb <- as_draws_matrix(fit_bdeb$fit$draws())
pars_compare <- c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L")

cat("\nBayesianDEB posterior:\n")
s_bdeb <- summarise_draws(draws_bdeb[, pars_compare],
  median = median,
  q5 = ~quantile(.x, 0.05),
  q95 = ~quantile(.x, 0.95),
  ess_bulk = ess_bulk)
print(as.data.frame(s_bdeb), digits = 3)
cat(sprintf("Wall-clock time: %.1f seconds\n", t_bdeb["elapsed"]))


## *******************************************************************
## 2. deBInfer
## *******************************************************************

cat("\n*** deBInfer fit ***\n")

## --- ODE function (must be written manually for deBInfer) ---

deb_ode <- function(t, state, pars) {
  E <- state["E"]
  V <- state["V"]
  L <- V^(1/3)

  p_Am  <- pars["p_Am"]
  p_M   <- pars["p_M"]
  kappa <- pars["kappa"]
  v     <- pars["v"]
  E_G   <- pars["E_G"]
  f     <- 1.0

  p_A <- f * p_Am * L^2
  p_C <- E * v * L / (E + E_G * V)
  p_Mf <- p_M * V

  dE <- p_A - p_C
  dV <- max((kappa * p_C - p_Mf) / E_G, 0)

  list(c(dE = dE, dV = dV))
}

## --- Observation model (manual likelihood) ---

deb_obs <- function(data, sim.data, samp) {
  sigma_L <- samp[["sigma_L"]]
  t_obs <- data$time
  L_obs <- data$length

  sol <- sim.data
  L_pred <- approx(sol[, "time"], sol[, "V"]^(1/3),
                   xout = t_obs)$y

  if (any(is.na(L_pred)) || any(L_pred <= 0)) return(-Inf)

  sum(dnorm(L_obs, L_pred, sigma_L, log = TRUE))
}

## --- Parameter specification for deBInfer ---

p_Am_de  <- debinfer_par(name = "p_Am", var.type = "de",
  fixed = FALSE, value = 5.0, prior = "lnorm",
  hypers = list(meanlog = 1.5, sdlog = 0.5),
  prop.var = 0.1, samp.type = "rw")

p_M_de   <- debinfer_par(name = "p_M", var.type = "de",
  fixed = FALSE, value = 0.5, prior = "lnorm",
  hypers = list(meanlog = -1.0, sdlog = 0.5),
  prop.var = 0.05, samp.type = "rw")

kappa_de <- debinfer_par(name = "kappa", var.type = "de",
  fixed = FALSE, value = 0.75, prior = "beta",
  hypers = list(shape1 = 3, shape2 = 2),
  prop.var = 0.02, samp.type = "rw")

v_de     <- debinfer_par(name = "v", var.type = "de",
  fixed = FALSE, value = 0.2, prior = "lnorm",
  hypers = list(meanlog = -1.5, sdlog = 0.5),
  prop.var = 0.05, samp.type = "rw")

E_G_de   <- debinfer_par(name = "E_G", var.type = "de",
  fixed = FALSE, value = 400, prior = "lnorm",
  hypers = list(meanlog = 6.0, sdlog = 0.5),
  prop.var = 20, samp.type = "rw")

sigma_de <- debinfer_par(name = "sigma_L", var.type = "obs",
  fixed = FALSE, value = 0.02, prior = "norm",
  hypers = list(mean = 0, sd = 0.05),
  prop.var = 0.005, samp.type = "rw")

# Init-state parameters must be named to match the ODE state vector
# (deb_ode reads state["E"], state["V"]).
E_de     <- debinfer_par(name = "E", var.type = "init",
  fixed = TRUE, value = 1.0)

V_de     <- debinfer_par(name = "V", var.type = "init",
  fixed = TRUE, value = 0.1^3)

## --- Prepare data for deBInfer ---

obs_data <- data.frame(time = df1$time, length = df1$length)

# de_mcmc() requires a debinfer_parlist (from setup_debinfer), not a
# plain list.
mcmc_pars <- setup_debinfer(p_Am_de, p_M_de, kappa_de,
                            v_de, E_G_de, sigma_de, E_de, V_de)

## --- Run deBInfer MCMC ---

t_debinfer <- system.time({
  mcmc_de <- tryCatch({
    de_mcmc(N = 20000, data = obs_data,
      de.model = deb_ode,
      obs.model = deb_obs,
      all.params = mcmc_pars,
      Tmax = max(df1$time),
      data.times = df1$time,
      cnt = 500,
      plot = FALSE,
      verbose.mcmc = FALSE)
  }, error = function(e) {
    cat("deBInfer fitting failed:", e$message, "\n")
    NULL
  })
})

if (!is.null(mcmc_de)) {
  ## Extract posterior (discard first 50% as burnin)
  burnin <- 10000
  chain <- as.data.frame(mcmc_de$samples[(burnin + 1):20000, ])

  cat("\ndeBInfer posterior:\n")
  for (p in pars_compare) {
    if (p %in% names(chain)) {
      vals <- chain[[p]]
      cat(sprintf("  %-10s: median = %.3f, 90%% CI = [%.3f, %.3f]\n",
        p, median(vals),
        quantile(vals, 0.05), quantile(vals, 0.95)))
    }
  }
  cat(sprintf("Wall-clock time: %.1f seconds\n",
              t_debinfer["elapsed"]))

  ## Effective sample size (using coda)
  if (requireNamespace("coda", quietly = TRUE)) {
    mcmc_obj <- coda::mcmc(as.matrix(chain[, pars_compare]))
    ess <- coda::effectiveSize(mcmc_obj)
    cat("\nEffective sample sizes:\n")
    print(round(ess))
    cat(sprintf("ESS/min (median across params): %.0f\n",
      median(ess) / (t_debinfer["elapsed"] / 60)))
  }
} else {
  cat("deBInfer comparison skipped (fitting failed)\n")
}


## *******************************************************************
## 3. Summary
## *******************************************************************

cat("\n*** Comparison summary ***\n")
cat(sprintf("BayesianDEB: %.0f sec, %d lines of model code\n",
            t_bdeb["elapsed"], 6))
cat(sprintf("deBInfer:    %.0f sec, ~80 lines of model code\n",
            t_debinfer["elapsed"]))

## BayesianDEB ESS
ess_bdeb <- summarise_draws(draws_bdeb[, pars_compare],
  ess_bulk = ess_bulk)$ess_bulk
cat(sprintf("\nBayesianDEB median ESS/min: %.0f\n",
  median(ess_bdeb) / (t_bdeb["elapsed"] / 60)))
