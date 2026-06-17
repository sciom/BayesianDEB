# BayesianDEB

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19443804.svg)](https://doi.org/10.5281/zenodo.19443804)
<!-- badges: end -->

**BayesianDEB** provides a Bayesian framework for Dynamic Energy Budget
(DEB) modelling in R, using [Stan](https://mc-stan.org/) as the
computational backend via
[cmdstanr](https://mc-stan.org/cmdstanr/).

The package implements the standard DEB model (Kooijman, 2010) as a
Bayesian state-space model.  Parameters are estimated via Hamiltonian
Monte Carlo (NUTS; Hoffman & Gelman, 2014), with full posterior
uncertainty propagation to derived biological quantities such as ultimate
length, von Bertalanffy growth rate, EC50, and NEC.

## Overview

Four model types cover the most common DEB applications:

| Type | States | Stan model | Use case |
|------|--------|------------|----------|
| `"individual"` | E, V | `bdeb_individual_growth` | Single organism growth |
| `"growth_repro"` | E, V, R | `bdeb_growth_repro` | Growth + reproduction |
| `"hierarchical"` | E, V + RE | `bdeb_hierarchical_growth` | Multi-individual, partial pooling |
| `"debtox"` | E, V, R, D | `bdeb_debtox` | Ecotoxicology (TKTD) |

Each model is bundled as a pre-written Stan program; the R layer handles
data preparation, prior specification, fitting, diagnostics, and
visualisation.

> **Note on model maturity.** The `"individual"` and `"growth_repro"`
> models have the most complete downstream tooling (L\_hat trajectories,
> posterior predictive checks, log\_lik for LOO-CV).  The
> `"hierarchical"` and `"debtox"` models use `reduce_sum` for
> within-chain parallelism, which means PPC and L\_hat are not stored in
> Stan; trajectory plots are instead produced via R-side forward
> simulation from the posterior.  All four model types support full
> diagnostics, summary, derived quantities, and trace/posterior/pairs
> plots.

## Feature status

| Feature | Status | Notes |
|---------|:------:|-------|
| Individual growth model | **stable** | 2-state DEB, full PPC |
| Growth + reproduction model | **stable** | 3-state DEB, NegBin/Poisson offspring |
| Hierarchical growth model | **stable** | Non-centred param., `reduce_sum` threading |
| DEBtox (TKTD) model | **beta** | Group-level (1 ODE/conc.); hierarchical DEBtox planned |
| Observation models (normal, lognormal, student-t) | **stable** | Switchable per model via `obs_*()` |
| Arrhenius temperature correction | **stable** | Global correction on rate parameters |
| Prior specification | **stable** | 6 families, weakly informative defaults |
| Species-specific priors (`prior_species()`) | **stable** | Built-in for *E. fetida*, *E. andrei*, *F. candida*, *D. magna*, *L. rubellus* |
| Derived quantities (L_m, L_inf, k_M, g, r_B) | **stable** | Kooijman (2010) formulas, dimensional tests |
| EC50 / NEC extraction | **stable** | Analytical, full posterior |
| Within-chain parallelism (`reduce_sum`) | **stable** | Hierarchical + DEBtox models |
| Trajectory plots for hierarchical/debtox | **stable** | R-side LSODA simulation from posterior |
| PPC for hierarchical/debtox | not available | Use trajectory plots and `bdeb_diagnose()` |
| Survival endpoint | planned | Not in current API; will be added in future version |
| Per-observation temperature (time-varying T) | planned | Currently global only |
| Multiple stress modes (maintenance, growth cost) | planned | Assimilation stress implemented |
| Prior-posterior comparison plot | **stable** | `plot(fit, type = "prior_posterior")` |
| Prior predictive checks | **stable** | `bdeb_prior_predictive()` with plot method |
| Standalone simulators | **stable** | `deb_simulate()`, `debtox_simulate()` |
| Maturity dynamics (wider) | planned | Currently post-puberty only |
| Hierarchical DEBtox (individual-level TKTD) | planned | Currently group-level only |
| Weight endpoint | planned | Requires shape coefficient in model |
| Shape coefficient ($\delta_M$) estimation | planned | Currently user-supplied conversion |
| Real-data benchmark dataset | planned | Simulated data validated against AmP |

## Current limitations

- `"individual"` and `"growth_repro"` support **exactly one individual**;
  multi-individual data require `"hierarchical"`.
- `"debtox"` fits **one ODE per concentration group**, not per individual.
  Multi-individual data within a group are aggregated to means (with warning).
  Hierarchical individual-level DEBtox is planned.
- `"hierarchical"` models random effects on `p_Am` only; other parameters
  are shared across individuals.
- Survival endpoint is not implemented (not in API).
- Temperature correction is global (single T for the entire experiment);
  time-varying temperature is planned.
- Food level (`f_food`) is constant; time-varying food is not supported.
- All lengths are structural ($L = V^{1/3}$); the shape coefficient
  $\delta_M$ is not estimated (user must convert physical lengths).
- R-side simulators (`deb_simulate`, trajectory plots, prior predictive)
  use the LSODA solver from `deSolve`, matching Stan's BDF solver.

## Installation

`BayesianDEB` installs from CRAN like any other package:

```r
install.packages("BayesianDEB")
```

Fitting models additionally requires **cmdstanr** and a working
**CmdStan** toolchain.  `cmdstanr` is not on CRAN, so install it from the
Stan r-universe and build CmdStan once:

```r
# 1. Install cmdstanr from the Stan r-universe (not on CRAN)
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# 2. Build CmdStan (one-time, ~10 min)
cmdstanr::install_cmdstan()
```

The package loads, prints and runs all non-fitting functions without
cmdstanr; `bdeb_fit()` checks for the toolchain at runtime and gives an
informative error if it is missing.  For the development version:

```r
# remotes::install_github("sciom/BayesianDEB")
```

## Quick start

```r
library(BayesianDEB)

# Bundled example data: 21 Eisenia fetida individuals, 12 weeks
data(eisenia_growth)

# 1. Prepare data (single individual)
dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])

# 2. Specify model with custom priors
mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am  = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa = prior_beta(a = 3, b = 2)
  ))

# 3. Fit
fit <- bdeb_fit(mod, chains = 4, iter_sampling = 1000)

# 4. Inspect
bdeb_diagnose(fit)
plot(fit, type = "trajectory")

# 5. Posterior predictive check
ppc <- bdeb_ppc(fit)
plot(ppc)

# 6. Derived quantities
bdeb_derived(fit, quantities = c("L_inf", "growth_rate"))
```

## Citation

If you use BayesianDEB in your work, please cite:

> Hackenberger, B.K., Djerdj, T. and Hackenberger, D.K. (2026).
> BayesianDEB: Bayesian Dynamic Energy Budget Modelling. R package
> version 0.2.1. doi:10.5281/zenodo.19443804.
> https://github.com/sciom/BayesianDEB

## Key references

- Kooijman, S.A.L.M. (2010). *Dynamic Energy Budget Theory for Metabolic
  Organisation*. 3rd ed. Cambridge University Press.
  doi:[10.1017/CBO9780511805400](https://doi.org/10.1017/CBO9780511805400)
- Jager, T., Heugens, E.H.W. and Kooijman, S.A.L.M. (2006). Making sense
  of ecotoxicological test results: towards application of process-based
  models. *Ecotoxicology*, 15(3), 305--314.
  doi:[10.1007/s10646-006-0060-x](https://doi.org/10.1007/s10646-006-0060-x)
- Carpenter, B. et al. (2017). Stan: A probabilistic programming language.
  *Journal of Statistical Software*, 76(1).
  doi:[10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

## License

MIT
