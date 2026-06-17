# BayesianDEB 0.2.1

Second-round JSS revision: diagnostics display polish, fully runnable
(toolchain-gated) examples, and replication / manuscript reproducibility
fixes.

* `print.bdeb_diagnostics()` now prints a compact table by default,
  hiding the per-time-point latent states (`x_sol[i,j]`, `L_hat[i]`);
  the new `full = TRUE` argument restores the complete table (also
  available via `summary(x)$table`).
* `plot.bdeb_diagnostics()` likewise hides the per-time-point latent
  states by default for a shorter, more readable R-hat / ESS plot;
  the new `full = TRUE` argument plots every monitored quantity. The
  subtitle notes how many latent rows were hidden.
* `bdeb_diagnose()` and the `bdeb_diagnostics` print/summary/plot
  methods, plus `bdeb_loo()`, `bdeb_derived()` and `bdeb_summary()`,
  moved their examples from `\dontrun{}` to a cmdstanr-gated
  `\donttest{}` block, so `example()` runs them when a CmdStan
  toolchain is available (mirroring `bdeb_fit()`).
* `plot_dose_response()` now draws posterior-median EC50 (dashed red)
  and NEC (dotted green) reference lines, drops non-finite draws, and
  clips the view with `coord_cartesian()` so degenerate draws no longer
  appear as vertical artefacts.
* Replication material: a single-command `replicate_all.R` reproduces
  every manuscript figure, table and printed number in ~1 min from
  cached draws; every manuscript code listing is now executed and
  printed; outputs use a fixed seed for bit-identical reproduction; and
  the README maps every figure to the script that produces it.

# BayesianDEB 0.2.0

JSS revision release.  Comprehensive S3 class system overhaul,
expanded methods coverage, and reproducible replication package.

## Breaking changes
* `bdeb_diagnose()` now returns a `bdeb_diagnostics` S3 object.
  All output goes through `print()`, `summary()`, and `plot()`
  (`type = "rhat"` or `"ess"`).  Direct list access still works
  (`n_divergent`, `summary`, ...).
* `bdeb_summary()` is deprecated; use `summary(fit)` on a
  `bdeb_fit` object instead.  The wrapper still works but emits
  a deprecation warning and will be removed in a future release.

## New methods
* `bdeb_data`: `print()`, `summary()`, `plot()`.
* `bdeb_model`: `print()`, `summary()`, `plot()`.
* `bdeb_prior`: `print()`, `summary()`, `plot()`.
* `bdeb_prediction`: `print()`, `summary()` (the latter returns a
  tidy `time / lower / median / upper` data frame).
* `bdeb_diagnostics`: `print()`, `summary()`, `plot()`.
* `bdeb_fit` (additional `lm`-style methods): `summary()` (now the
  primary posterior-summary API, accepts `pars` and `prob`),
  `confint()` (posterior credible intervals), `fitted()` (posterior
  median/mean of \eqn{\hat{L}_i}), `residuals()` (observed minus
  fitted), `nobs()` (observation count), `vcov()` (posterior
  covariance of model parameters), and `logLik()` (log-pointwise
  predictive density, lppd).  `fitted()`, `residuals()`, and
  `logLik()` are available for `"individual"` and `"growth_repro"`
  models only.
* `bdeb_derived()` is now an S3 generic with a `bdeb_fit` method
  (`bdeb_derived.bdeb_fit`).  Existing calls are unchanged; the
  dispatch enables future support for derived quantities on
  prior-only or simulated objects.

## Bug fixes
* `bdeb_diagnose()` no longer fails on real-data fits whose
  generated quantities (`log_lik`, `L_rep`) contain NaN draws
  from sporadic ODE-solver failures.  The `summarise_draws()`
  quantile lambdas now pass `na.rm = TRUE`.
* `plot()` on a `bdeb_fit` (`type = "trace"`, `"posterior"`,
  `"pairs"`) now subsets draws to the requested parameters
  before delegating to `bayesplot`.  Without this, NaN draws in
  unrelated generated quantities caused
  `bayesplot::prepare_mcmc_array()` to abort with
  `"NAs not allowed in 'x'"`.

## Stan
* Increased `ode_bdf_tol` `max_num_steps` from 1e4 to 1e5 in all
  four Stan models to reduce CVode mxstep messages during warmup.

## Vignettes
* Conditional execution via `requireNamespace("cmdstanr")`.
* Fixed `plot(fit, type = "pairs")` returning `NULL`.

## Replication material
* Reorganised JSS replication material into a single zip with
  `README`, `data/`, `outputs/`, and `lite`/`full` execution modes.
* Bundled `curves.txt` locally to remove external dependency.

---

# BayesianDEB 0.1.4

CRAN-review compliance release.  Addresses reviewer feedback on the use
of `T` as an identifier, which can shadow R's built-in `T` symbol
(= `TRUE`).

## Breaking changes
* `arrhenius()`: first argument renamed from `T` to `temp`.  Positional
  calls are unaffected (`arrhenius(298.15)`).  Code that passed the
  argument by name (`arrhenius(T = 298.15)`) must be updated to
  `arrhenius(temp = 298.15)`.
* `bdeb_model(temperature = ...)`: the list field `T` was renamed to
  `T_obs` to match the Stan data naming and to remove the `T` shadow.
  Replace `list(T = 298.15, T_ref = ..., T_A = ...)` with
  `list(T_obs = 298.15, T_ref = ..., T_A = ...)`.

## CRAN compliance
* `bdeb_diagnose()` and `bdeb_ec50()` now expose a `verbose = TRUE`
  argument.  All user-facing output (diagnostic alerts and summary
  tables) is routed through `cli` / [message()] rather than direct
  `print()` calls, so it can be silenced with [suppressMessages()] or
  by passing `verbose = FALSE`.  Return values are unchanged.

## Documentation
* Updated `R/utils.R`, `R/data_prep.R`, `R/model_spec.R` and
  `man/arrhenius.Rd`, `man/bdeb_model.Rd`,
  `man/temperature_to_stan_data.Rd`,
  `man/build_stan_data_individual.Rd` to reflect the renaming.
* Case-study vignette updated to use `arrhenius(temp = ...)`.
* Tests updated for the new `temperature$T_obs` field.
* `\dontrun{}` replaced with runnable examples for `prior_species()`,
  `bdeb_tox()` and `bdeb_prior_predictive()` (all execute in < 0.2 s
  against bundled datasets).  `bdeb_fit()` retains `\dontrun{}`
  because it requires the external CmdStan toolchain and a single
  Stan compilation + MCMC run takes > 30 seconds; the Rd comment
  explains the reason.

## Citation
* Switched README badge and `CITATION.cff` from the version-specific
  DOI (`10.5281/zenodo.19500753`, v0.1.3) to the **concept DOI**
  (`10.5281/zenodo.19443804`), which always resolves to the latest
  archived version.  The v0.1.4 version DOI is
  `10.5281/zenodo.19642839` for anyone needing to cite this specific
  release.

---

# BayesianDEB 0.1.3

New features and data release.

## New features
* `prior_species()`: species-specific priors from the AmP collection for
  *E. fetida*, *E. andrei*, *F. candida*, *D. magna*, and *L. rubellus*.
* `plot(fit, type = "prior_posterior")`: prior vs. posterior density
  comparison plot.

## New datasets
* `eisenia_neuhauser`: real *E. fetida* growth data from Neuhauser,
  Hartenstein & Kaplan (1980), 37 group-mean measurements over 250 days.
* `eisenia_cd`: real *Eisenia andrei* cadmium toxicity data from
  Van Gestel et al. (1991), 5 concentration groups over 85 days.

## Documentation
* Updated DOI to Zenodo 10.5281/zenodo.19500753.

---

# BayesianDEB 0.1.2

ODE solver upgrade and bug fix release.

## Breaking changes
* R-side ODE integration now uses `deSolve::lsoda()` (adaptive
  BDF/Adams) instead of fixed-step Euler. This matches the BDF
  solver used in Stan, ensuring numerical consistency between
  R-side simulation and Stan-side inference. The `dt` parameter
  now controls output resolution, not integration accuracy.
* `deSolve` (>= 1.40) added to Imports.

## Bug fixes
* Fixed `bdeb_ec50()$NEC` returning a data frame instead of numeric
  draws — `median(ec50$NEC)` now works correctly.
* Fixed `build_stan_data_hierarchical()` crashing when growth data
  starts at `time = 0` (same fix as individual model in 0.1.1).
* Fixed `build_stan_data_debtox()` crashing when growth or
  reproduction data starts at `time = 0`.
* Fixed `build_stan_data_growth_repro()` t=0 handling producing NA
  index matches — replaced convoluted fix with consistent epsilon
  shift across growth and reproduction times.

## Internal
* Renamed internal helpers `sim_deb_euler()` / `sim_debtox_euler()`
  to `sim_deb_lsoda()` / `sim_debtox_lsoda()` to reflect solver.
* Updated package documentation to describe LSODA solver instead
  of Euler.

---

# BayesianDEB 0.1.1

Bug fix and hardening release.

## Bug fixes
* Fixed `ode_bdf_tol` crash when growth data starts at `time = 0`:
  replaced with `1e-3` epsilon since Stan requires `t_obs > t0`.
* Fixed `%||%` operator not imported: replaced with explicit
  `if (is.null(...))` in `bdeb_predict()`.
* Fixed silent `rnorm()` fallback for unknown prior families in
  `bdeb_prior_predictive()`: now throws informative error.
* Fixed `bdeb_predict(newdata = NULL)` using `t_obs` instead of
  `t_L` for `"growth_repro"` models.
* Fixed `bdeb_predict()` silently returning raw draws when `L_hat`
  is missing: now throws informative error.
* Fixed DEBtox reproduction using cumulative `R` instead of interval
  `delta_R = k_R * (R(t_end) - R(t_start))`: now consistent with
  `growth_repro` model.
* Fixed `phi_R` prior hardcoded in DEBtox Stan model: now uses
  `prior_phi_R_mu/sd` from R prior system.
* Fixed `deb_fluxes()` computing `p_G` with redundant `/ E_G * E_G`.
* Fixed `growth_rate` formula: was `v/3 * p_M/(kappa*E_G)`, now
  correct Kooijman Eq. 3.23: `k_M * g / (3*(f+g))`.

## Safety improvements
* `"individual"` and `"growth_repro"` models now hard-error (not warn)
  when data contains multiple individuals.
* DEBtox auto-aggregates multi-individual concentration groups to
  means with explicit warning.
* Removed `survival` argument from `bdeb_data()` (was accepted but
  not implemented — "feature mirage").
* All prior constructors validate hyperparameters (sigma > 0, a > 0,
  rate > 0, nu > 1).
* `bdeb_data()` validates `f_food` in [0, 1].
* `bdeb_fit()` validates all sampling parameters (chains, iterations,
  adapt_delta, max_treedepth, threads).
* `arrhenius()` validates T > 0, T_ref > 0, T_A >= 0.
* `bdeb_model()` validates prior objects, observation objects, and
  temperature fields (must be positive finite scalars).
* `bdeb_fit()` wraps Stan compilation and sampling in `tryCatch`
  with contextualised error messages.
* `repro_to_intervals()` warns when dropping individuals with < 2
  observations.
* `validate_growth()` warns when max length > 10 cm (possibly
  physical rather than structural length).
* DEBtox reproduction time matching now uses strict `match()` instead
  of nearest-neighbour `which.min()`.

## New features
* `bdeb_loo()`: LOO cross-validation via `loo::loo()` with
  `endpoint` argument for `"growth_repro"` models.
* `bdeb_prior_predictive()`: R-side prior predictive simulation
  with `print()` and `plot()` methods.
* `bdeb_session_info()`: reproducibility report (R, Stan, package
  versions, fit configuration, Stan model hash).
* `coef.bdeb_fit()`: S3 method returning posterior medians or means.
* `deb_simulate()`, `debtox_simulate()`: standalone DEB/DEBtox
  simulators, independent of Stan.
* `L_m` added to `bdeb_derived()`: maximum structural length at f=1.
* `plot.bdeb_prediction()`: S3 plot method for prediction objects.
* Observation model switching now fully implemented in all 4 Stan
  models via integer flags (no recompilation needed).
* Arrhenius temperature correction implemented in all 4 Stan models.
* Within-chain parallelism via `reduce_sum` for hierarchical and
  DEBtox models.
* `threads_per_chain` argument in `bdeb_fit()`.
* `seed` argument in `bdeb_predict()`, `plot_dose_response()`, and
  `plot(fit, type = "trajectory")`.
* `dt` and `t_end` arguments in `plot_dose_response()`.
* Reproducibility metadata (seed, versions, timestamp) stored in
  `bdeb_fit` object.
* `inst/CITATION` with package and Kooijman (2010) entries.
* `CITATION.cff` with Zenodo DOI.

## Documentation
* Structural vs physical length documented throughout (delta_M).
* All derived quantity formulas reference Kooijman (2010) equations.
* Prior calibration rationale documented with AmP ranges.
* Feature status table in README with stable/beta/planned.
* "Current limitations" section in README.
* "Numerical layers" section in package docs (Stan exact vs R-side
  approximate).
* Lifecycle annotations on all functions and models.

## Tests
* 901 tests across 16 files (up from 58).
* New test categories: API contracts, scientific consistency,
  snapshot/regression, deep validation, end-to-end integration.

---

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
* Sensible weakly informative default priors via `prior_default()`,
  calibrated against the AmP collection (Marques et al., 2018).
* Observation model selection: `obs_normal()`, `obs_lognormal()`,
  `obs_student_t()`, `obs_poisson()`, `obs_negbinom()`.
* MCMC fitting via `cmdstanr` with `bdeb_fit()`.
* Convergence diagnostics: `bdeb_diagnose()` reports R-hat, ESS,
  divergences, and E-BFMI (Vehtari et al., 2021).
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

## Vignettes
* `getting_started`: overview of the package workflow.
* `case_study_eisenia_folsomia`: full workflow with Eisenia growth
  (individual + hierarchical) and Folsomia DEBtox analysis.
