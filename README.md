# BayesianDEB

<!-- badges: start -->
<!-- badges: end -->

**BayesianDEB** provides a Bayesian framework for Dynamic Energy Budget
(DEB) modelling in R, using [Stan](https://mc-stan.org/) as the
computational backend via
[cmdstanr](https://mc-stan.org/cmdstanr/).

## Overview

The package implements four model types that cover common DEB applications:
individual growth, growth with reproduction, hierarchical multi-individual
analysis, and ecotoxicological (DEBtox/TKTD) models. Each model is
bundled as a pre-written Stan program; the R layer handles data
preparation, prior specification, fitting, diagnostics, and visualisation.

## Installation

BayesianDEB requires **cmdstanr** and a working **CmdStan** installation.

```r
# 1. Install cmdstanr from r-universe
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# 2. Install CmdStan (one-time)
cmdstanr::install_cmdstan()

# 3. Install BayesianDEB
# install.packages("BayesianDEB")
# or from source:
# remotes::install_github("bhackenberger/BayesianDEB")
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

## Model types

| Type | States | Stan model | Use case |
|------|--------|------------|----------|
| `"individual"` | E, V | `bdeb_individual_growth` | Single organism growth |
| `"growth_repro"` | E, V, R | `bdeb_growth_repro` | Growth + reproduction |
| `"hierarchical"` | E, V + RE | `bdeb_hierarchical_growth` | Multi-individual, partial pooling |
| `"debtox"` | E, V, R, D | `bdeb_debtox` | Ecotoxicology (TKTD) |

## Citation

If you use BayesianDEB in your work, please cite:

> Hackenberger, B.K. (2025). BayesianDEB: Bayesian Dynamic Energy Budget
> Modelling. R package version 0.1.0.

## License

GPL (>= 3)
