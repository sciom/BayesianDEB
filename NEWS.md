# BayesianDEB 0.1.0

Initial release.

## Models
* Individual-level growth model (2-state DEB: reserve, structure).
* Growth + reproduction model (3-state DEB: reserve, structure,
  reproduction buffer) with negative binomial observation model.
* Hierarchical multi-individual growth model with non-centred
  parameterisation and partial pooling of assimilation rates.
* DEBtox (TKTD) model for ecotoxicology with scaled internal damage,
  stress on assimilation, and analytical EC50/NEC computation.

## Features
* Declarative model specification via `bdeb_model()`.
* Prior specification functions for all standard DEB parameters:
  `prior_lognormal()`, `prior_beta()`, `prior_halfnormal()`, and others.
* Sensible weakly informative default priors via `prior_default()`.
* Observation model selection: `obs_normal()`, `obs_lognormal()`,
  `obs_student_t()`, `obs_poisson()`, `obs_negbinom()`.
* MCMC fitting via `cmdstanr` with `bdeb_fit()`.
* Convergence diagnostics: `bdeb_diagnose()` reports R-hat, ESS,
  divergences, and E-BFMI.
* Posterior predictive checks with `bdeb_ppc()`.
* Derived quantity computation: `bdeb_derived()` for ultimate length,
  von Bertalanffy growth rate, and somatic maintenance rate constant.
* Publication-quality plots: trace, posterior density, pairs, trajectory,
  and PPC overlays.
* DEBtox helpers: `bdeb_tox()`, `bdeb_ec50()`, `plot_dose_response()`.
* Utility functions: `arrhenius()` temperature correction, `deb_fluxes()`
  energy flux calculator, `repro_to_intervals()` data converter.

## Data
* `eisenia_growth`: simulated *Eisenia fetida* growth (21 individuals,
  12 weeks).
* `folsomia_repro`: simulated *Folsomia candida* reproduction test
  (5 Cd concentrations, 6 replicates).
* `debtox_growth`: simulated growth under toxicant exposure
  (4 concentrations, 10 individuals each, 6 weeks).
