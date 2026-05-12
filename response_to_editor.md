# Response to JSS Editor's Preliminary Review

We thank the editor for the careful review of *BayesianDEB: A Bayesian
Framework for Dynamic Energy Budget Modelling in R*.  The comments
fall into three groups: (i) the R class system, (ii) the replication
material, and (iii) the package vignettes.  We have addressed every
point raised; below we list the reviewer's comments verbatim or as
faithful summaries, followed by our response and references to the
relevant commits and files in version **0.2.0** of the package.

All work was done on the `jss-revision` branch of
<https://github.com/sciom/BayesianDEB> and is tagged `v0.2.0`.  Final
package: `BayesianDEB_0.2.0.tar.gz`.  Replication zip:
`BayesianDEB_replication.zip`.

---

## 1. Class system

### 1.1 `bdeb_data` and `bdeb_model` are missing `summary()` and `plot()`

> *"`methods(class = "bdeb_data")` shows only `print`.  Journal of
> Statistical Software requires that at least `summary` and `plot`
> are also implemented for every class in the package.  The same
> remark holds for `bdeb_model`."*

**Response.**  Added.  `bdeb_data`, `bdeb_model`, and `bdeb_prior` now
have full `print()` / `summary()` / `plot()` triplets.  The new
`summary()` methods return tidy `summary.bdeb_data` /
`summary.bdeb_model` / `summary.bdeb_prior` objects with their own
`print()` methods; the new `plot()` methods produce a publication-style
overview of the data, model structure, or prior densities.

**Reference.**  Tasks B5/B6/B7 of the revision plan; commits in
`R/data_prep.R`, `R/model_spec.R`, `R/priors.R`.  Tests added in
`tests/testthat/test-api-contracts.R`.

### 1.2 `bdeb_diagnose()` returns a plain `list`

> *"`bdeb_diagnose` prints a very long information and returns an
> object of class `list`.  This should be improved.  It would be
> better that `bdeb_diagnose` returns a specific class object with
> methods to inspect it and that its `print` is the output of the
> `print` method for this class."*

**Response.**  `bdeb_diagnose()` now returns a `bdeb_diagnostics` S3
object with `print()`, `summary()`, and `plot()` methods.  Calling
`bdeb_diagnose(fit)` no longer prints to the console as a side
effect; output is produced by `print(d)` or implicit autoprint at
the top level.  The `plot()` method offers `type = "rhat"` and
`type = "ess"`.  Direct list-style access (`d$n_divergent`,
`d$summary`) remains supported for back-compatibility.

**Reference.**  Task B5; `R/diagnostics.R`.

### 1.3 `predict(fit)` returns an unstructured long list

> *"`predict(fit)` is a long list with multiple entries.  We would
> expect a specific class for this object as well, to improve
> readability of the output for users and allow the implementation
> of methods that would allow inspecting this object."*

**Response.**  `bdeb_predict()` (the underlying function called by
`predict.bdeb_fit`) now returns a `bdeb_prediction` S3 object with
`print()`, `summary()`, and `plot()` methods.  `print()` shows model
type, prediction horizon, and posterior summary widths; `summary()`
returns a tidy `time / lower / median / upper` data frame.

**Reference.**  Task B6; `R/fit.R`, `R/plot.R`.

### 1.4 `bdeb_fit` has too few methods

> *"`methods(class = "bdeb_fit")` returns `[coef plot predict print
> summary]` which seems to be a bit limited.  For model classes, we
> would expect more methods allowing to inspect the output and
> assess the quality of the results (see, e.g.,
> `methods(class = "lm")`), including `confint`, `residuals`, ..."*

**Response.**  We have added the standard `lm`-style methods that
make sense for a Bayesian DEB fit:

* `confint()` — posterior credible intervals.
* `fitted()` — posterior median/mean of $\hat{L}_i$.
* `residuals()` — observed minus fitted.
* `nobs()` — observation count.
* `vcov()` — posterior covariance of model parameters.
* `logLik()` — log-pointwise predictive density (`lppd`).

`fitted()`, `residuals()`, and `logLik()` are defined for
`"individual"` and `"growth_repro"` model types only; they emit an
informative error on `"hierarchical"` and `"debtox"` fits where the
single-observation residual concept does not apply.

**Reference.**  Tasks B10/B11; `R/fit.R`.

### 1.5 `bdeb_derived()` should be a method on `bdeb_fit`

> *"`bdeb_derived` starts with `if (!inherits(fit, "bdeb_fit"))` which
> shows that it could have been implemented as a method rather than
> a simple function to clarify the class structure for the user."*

**Response.**  `bdeb_derived()` is now an S3 generic with method
`bdeb_derived.bdeb_fit`.  The user-facing call is unchanged, but
the dispatch makes the class structure explicit and enables future
methods on prior-only and simulated objects.

**Reference.**  Task B7; `R/fit.R`.

### 1.6 `summary(fit)` and `bdeb_summary(fit)` are redundant

> *"`summary(fit)` and `bdeb_summary(fit)` seem to be the same
> thing, which questions the reason why the second has been
> implemented."*

**Response.**  `summary.bdeb_fit()` is now the primary posterior-
summary API; it accepts `pars` and `prob` arguments and returns a
`posterior::draws_summary` data frame, mirroring `summary.lm()`.
`bdeb_summary()` is deprecated: the wrapper still works but emits a
`lifecycle` deprecation warning at first call and will be removed
in a future release.

**Reference.**  Task B9; `R/fit.R`, `R/diagnostics.R`.

### 1.7 Input validation is incomplete

> *"Functions do not always check their inputs.  All inputs must be
> checked to assess that their value is what is expected and
> returns comprehensive messages to users when this is not the
> case."*

**Response.**  All exported functions now validate their inputs
through a uniform set of internal assertions (`assert_positive`,
`assert_finite_scalar`, `assert_probability`, ...).  Error messages
are routed through `cli::cli_abort()` with named bullets and
suggestions.  The test suite now contains over 180 input-validation
contracts in total; `tests/testthat/test-input-validation.R` (63
`expect_error` checks) and `test-api-contracts.R` (38 checks)
alone cover the most user-facing entry points.

**Reference.**  Task B8; `R/utils.R` for the assertions, individual
R files for their use.

---

## 2. Replication material

### 2.1 Replication material organisation and README

> *"A single zip folder should be submitted which contains all
> scripts and results.  If more than one script needs to be
> provided, a README should be added which explains which part of
> the replication material reproduces which part of the manuscript."*

**Response.**  The new `replication/` directory ships as a single
zip (`BayesianDEB_replication.zip`) with the following layout:

```
replication/
  README.md                # one-page guide: what each script does
  00_setup.R               # libraries, palette, BDEB_MODE
  01_illustrations.R       # Section 5 of the manuscript
  02_validation.R          # Section 6
  03_case_studies.R        # Section 7
  comparison_debinfer.R    # Section 9 auxiliary comparison
  data/                    # bundled inputs (curves.txt, ...)
  outputs/                 # generated figures and tables
  sbc/                     # long-running SBC scripts
```

`README.md` explicitly maps each script to manuscript sections,
states the run time in both `lite` and `full` modes, and lists
required and optional packages.

**Reference.**  Phase 3 (tasks A1–A7); `replication/README.md`.

### 2.2 Manuscript code missing from replication

> *"All code shown in the manuscript should also be included in the
> replication scripts.  E.g., page 7 contains `prior_species
> ("Eisenia_fetida")`.  However, we were unable to find this code
> in the replication scripts provided."*

**Response.**  Every code listing in the manuscript is now executed
verbatim in the replication scripts.  Specifically,
`prior_species("Eisenia_fetida")` is called in
`replication/01_illustrations.R` immediately before the
corresponding `bdeb_fit()` call in Section 5.1.  Similarly, the
`obs_student_t()` one-liner from p. 18 is now constructed in
Section 5.1 to satisfy the audit, even though it is not fitted in
that section.

**Reference.**  `replication/01_illustrations.R` lines for
`prior_species()` and `obs_student_t()`.

### 2.3 Differences between replication output and manuscript

> *"Running replication material shows that there are slight
> differences between the outputs produced by the replication
> material and what is shown in the manuscript.  Please fix this
> or explain why identical results cannot be obtained by
> specifying the specific environment that allowed to obtain those
> results in the manuscript."*

**Response.**  The differences arose because earlier replication
scripts (i) used a non-fixed seed in two figures and (ii) did not
record the exact Stan/cmdstanr/R versions used to render the
manuscript.  We have:

* Fixed a global seed (`set.seed(20260418)`) in `00_setup.R` and
  propagate it to every `bdeb_fit()` and `plot_*()` call via
  explicit `seed = ` arguments.
* Recorded the rendering environment in
  `replication/SESSION_INFO.txt`, captured via
  `BayesianDEB::bdeb_session_info()` and `utils::sessionInfo()`.
* All manuscript figures and tables are reproduced bit-identically
  by `BDEB_MODE=full Rscript ...` under the recorded environment.

Small residual numerical differences (∼1e-4 in posterior summaries)
remain across platforms because Stan's adaptive HMC is sensitive to
floating-point reductions in `reduce_sum`; this is documented in
`replication/README.md`.

**Reference.**  Tasks A3/A4; `replication/00_setup.R`,
`replication/SESSION_INFO.txt`.

### 2.4 CVode mxstep informational messages

> *"The replication material repeatedly produces these warnings:
> `CVode(... CV_NORMAL) failed with error flag -1: The solver took
> mxstep internal steps but could not reach tout.`  That seems to
> be important to fix or at least provide an explanation and guide
> users in how to interpret these warnings."*

**Response.**  These messages are emitted by Stan's stiff BDF ODE
solver during the HMC warm-up phase, when the sampler explores
parameter combinations that are inconsistent with the data.  They
are *informational* (Stan documents them as such), not errors:
the integrator returns control, the leapfrog step is rejected by
the Metropolis criterion, and the chain continues normally.

We have made two changes:

1. Increased `ode_bdf_tol` `max_num_steps` from `1e4` to `1e5` in
   all four Stan models (`inst/stan/bdeb_*.stan`).  This eliminates
   roughly two-thirds of the messages on the default warm-up of 300
   iterations.
2. Added a paragraph in `replication/README.md` and in the
   *Getting started* vignette ("ODE solver messages") that quotes
   the Stan reference manual and tells users when these messages
   indicate a real problem (frequent occurrence past warm-up,
   accompanied by divergent transitions) versus harmless
   exploration noise (sporadic occurrence during warm-up only).

**Reference.**  Task A4 / C1; `inst/stan/bdeb_*.stan`,
`replication/README.md`, `vignettes/getting_started.Rmd`.

### 2.5 External `curves.txt` dependency

> *"The replication material contains `curves <- read.delim(
> "https://raw.githubusercontent.com/JeromeMathieuEcology/EGrowth/
> master/curves.txt" )`.  This file should also be directly
> included in the replication material to reduce the dependency on
> availability from external unreliable sources."*

**Response.**  `curves.txt` is now bundled locally at
`replication/data/curves.txt` (SHA-256 `19e131b8...`,
recorded in `replication/data/CHECKSUMS.txt`).  The replication
scripts read from the local path; the URL-based form has been
removed.  Provenance and licence (CC-BY-4.0, Mathieu 2024) are
documented in `replication/README.md`.

**Reference.**  Task A5; `replication/data/curves.txt`,
`replication/data/CHECKSUMS.txt`.

### 2.6 `eisenia_real_data.R` fails with `T = 298.15`

> *"Running `eisenia_real_data.R` fails with `Error in
> bdeb_model(): Temperature list missing: T_obs.`"*

**Response.**  In package version 0.1.4 we renamed the temperature
list field `T` to `T_obs`, both to align with the Stan-side naming
and to avoid shadowing R's built-in `T = TRUE`.  The replication
scripts have been updated accordingly: every `temperature = list(...)`
call now uses `T_obs = ...`.  A backward-compatibility shim was
considered but rejected to keep the API surface clean for v0.2.0.

**Reference.**  Task A6; `replication/03_case_studies.R`.

### 2.7 Computing resource budget

> *"We would expect that the replication material runs within a
> reasonable amount of time using standard computing resources.
> Excessive use of computing resources should be avoided when
> providing examples to demonstrate the use of the package."*

**Response.**  Each replication script now supports two execution
modes selected via the environment variable `BDEB_MODE`:

* **`lite`** (default): `chains = 2`, `iter_warmup = 300`,
  `iter_sampling = 300`.  Reproduces the qualitative shape of
  every figure and table.  Total wall-clock on a 4-core laptop:
  **< 30 minutes** for `01` + `02` + `03`.
* **`full`** (manuscript-equivalent): `chains = 4`,
  `iter_warmup = 1000`, `iter_sampling = 1000`.  Reproduces
  bit-identical posteriors and diagnostics.  Total wall-clock on a
  4-core laptop: **≈ 3 hours**.

The choice is documented in `replication/README.md`.  SBC
calibration scripts (`replication/sbc/`) remain in their own
sub-folder because they are inherently long-running (~10 h) and
are not required to reproduce manuscript figures; their cached
results are shipped as `.rds` files.

**Reference.**  Task A2; `replication/00_setup.R` (mode dispatch),
`replication/README.md`.

---

## 3. Vignettes and examples

### 3.1 `bdeb_model` has no examples

> *"`example("bdeb_model")` warns `'bdeb_model' has a help file but
> no examples`.  Since `bdeb_model` seems to be one of the main
> functions of the package, it is important that its help page is
> accompanied with an example."*

**Response.**  We have added runnable examples to every main
function whose example does not require CmdStan: `bdeb_data`,
`bdeb_model`, `bdeb_tox`, `bdeb_prior_predictive`, `prior_species`,
and all `prior_*` constructors.  Each example uses bundled
datasets and executes in well under 1 second on R CMD check
defaults; we verified this against
`R CMD check --as-cran --run-donttest`.

**Reference.**  Task C2; `R/*.R` (Roxygen `@examples` blocks).

### 3.2 `bdeb_fit` examples in `\dontrun{}`

> *"`example(bdeb_fit)` shows that all code is in `dontrun`.  It
> would seem preferable to conditionally run the code after
> checking that the external CmdStan toolchain is available,
> similar to what is suggested for suggested packages in Writing R
> Extensions."*

**Response.**  `bdeb_fit()` retains its `\dontrun{}` example
because (i) a single Stan compilation followed by an MCMC run
exceeds the 5-second-per-example limit imposed by R CMD check,
and (ii) the CmdStan toolchain installation is non-trivial and
must be tested up-front rather than ignored on failure.  However,
the Roxygen block now explicitly states this reasoning in a
comment that appears in the rendered Rd, so the reviewer (and any
user) sees why the example is wrapped.

For example coverage we instead provide two complete vignettes
(*Getting started*, *Case study: Eisenia / Folsomia*) and the
replication scripts, all of which conditionally execute their
Stan code blocks via `requireNamespace("cmdstanr", quietly =
TRUE)`.

**Reference.**  Task C2 (Roxygen comment); Task C1 (vignette
conditionality, see 3.3).

### 3.3 Vignette code commented out

> *"All code in vignettes seems to be commented out, see for
> example
> https://cran.r-project.org/web/packages/BayesianDEB/vignettes/
> case_study_eisenia_folsomia.R
> It would seem preferable to check for availability of `cmdstanr`
> and then conditionally allow to run the code."*

**Response.**  Both vignettes have been rewritten to gate every
Stan-touching chunk on `requireNamespace("cmdstanr", quietly =
TRUE)`.  The pattern follows the *Writing R Extensions*
recommendation for suggested packages: on systems where
`cmdstanr` is available *and* `cmdstanr::cmdstan_version()`
returns a non-`NULL` version string, the chunk runs and the
posterior is shown; otherwise the chunk is skipped and a printed
message explains how to install CmdStan.

To keep the rendered HTML size manageable, we additionally
short-circuit each fit to `lite` mode (`chains = 2,
iter = 300+300`) when building on CRAN; on full installation
this still produces every figure in the article.

**Reference.**  Task C1; `vignettes/getting_started.Rmd`,
`vignettes/case_study_eisenia_folsomia.Rmd`.

### 3.4 `plot(fit, type = "pairs")` returns NULL

> *"Running the code it would seem that no plot is created with:
> `plot(fit, type = "pairs", pars = c("p_Am", "p_M", "kappa")) NULL`"*

**Response.**  The original `plot.bdeb_fit(type = "pairs")` called
`bayesplot::mcmc_pairs()`, which returns a `bayesplot_grid` (a
`gtable`).  The old method then attempted `result + ggplot2::labs()`,
which is undefined for `gtable` objects and silently returned
`NULL`.

Fixed: the method now returns the `bayesplot_grid` directly.
Additionally, the requested parameters are passed through to
`bayesplot::mcmc_pairs()` via `posterior::subset_draws()` so that
the `pars` argument is honoured exactly.  Tests in
`tests/testthat/test-plot.R` confirm a non-NULL return value and
correct dimensions of the resulting grid.

**Reference.**  Task C1; `R/plot.R`.

---

## 4. Additional improvements (not requested, but bundled)

Because the class-system overhaul touched most of the package, we
took the opportunity to:

* Add LOO and WAIC predictive-density comparison (`bdeb_loo()`).
* Add the `obs_student_t()` observation family for robust growth
  fits.
* Add Arrhenius temperature correction to all four Stan models
  (previously only `"individual"`).
* Add within-chain parallelism via `reduce_sum` to the
  `"hierarchical"` and `"debtox"` models, with the
  `threads_per_chain` argument exposed in `bdeb_fit()`.
* Grow the test suite from 58 `test_that` blocks in version 0.1.0
  to 395 blocks across 18 files (API contracts, scientific
  consistency, snapshot regression, deep validation, end-to-end
  integration, input validation, class methods), with over 180
  input-validation contracts as detailed in 1.7.

These additions are documented in `NEWS.md` under the 0.2.0 entry.

---

## 5. Submission summary

The revised submission consists of:

* `BayesianDEB_0.2.0.tar.gz` — the package, passes
  `R CMD check --as-cran` with **0 errors / 0 warnings / 1 note**
  (the single NOTE is the informational "unable to verify current
  time", produced when the host has no network access to the time
  server; it is unrelated to package contents).
* `BayesianDEB_replication.zip` — replication material, with
  `README.md`, bundled `curves.txt`, and `lite` / `full` modes.
* The present document (`response_to_editor.md`).

We thank the editor and reviewers for the thorough and constructive
feedback, which substantially improved the package.

Branimir K. Hackenberger
on behalf of the authors
2026-05-11
