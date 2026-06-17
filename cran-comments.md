## R CMD check results

0 errors | 0 warnings | 1 note

Local: Ubuntu 22.04.5 LTS, R 4.5.2, CmdStan 2.36.0.

## Submission of 0.2.1

This release brings the CRAN version in line with the version submitted
to the Journal of Statistical Software (JSS requires that the published
package and the manuscript's package version match; CRAN currently holds
0.1.4).  Version 0.2.0 added a full S3 class system (print/summary/plot
for every exported class), additional `lm`-style methods on `bdeb_fit`
(`confint`, `fitted`, `residuals`, `nobs`, `vcov`, `logLik`), LOO/WAIC
comparison, and an expanded test suite.  Version 0.2.1 (this submission)
adds the second-round JSS revision: all exported examples are now
cmdstanr-gated `\donttest{}` blocks (no `\dontrun{}` remains), the
`bdeb_diagnostics` print and plot methods hide per-time-point latent
states by default for a compact display, and the replication material
ships a single `replicate_all.R` reproducing every manuscript result in
about one minute.  See `NEWS.md` for details.

## History (0.1.3 -> 0.1.4)

The following issues were fixed in 0.1.4:

1. **No identifier named `T` or `F`.** The argument `T` in `arrhenius()`
   was renamed to `temp` to avoid shadowing R's built-in `T` symbol
   (= `TRUE`); the list field `temperature$T` in `bdeb_model()` was
   renamed to `T_obs` for the same reason and to match the Stan-side
   naming.  No `T`/`F` is used as a logical literal anywhere in the
   package.

2. **Runnable examples instead of `\dontrun{}`.**  Examples in
   `prior_species()`, `bdeb_tox()` and `bdeb_prior_predictive()` were
   rewritten to use bundled datasets and now execute in < 0.2 s each.
   `bdeb_fit()` retains `\dontrun{}` because it requires the external
   CmdStan toolchain (not on CRAN) and a single fit takes > 30 seconds
   (Stan compilation + MCMC).  An explanatory comment precedes the
   `\dontrun{}` block in the Rd file.

3. **Suppressible output.**  Direct `print()` calls inside
   `bdeb_diagnose()` and `bdeb_ec50()` were removed; all user-facing
   output now goes through `cli` / `message()` and can be silenced
   with `suppressMessages()`.  Both functions also gained a new
   `verbose = TRUE` argument for explicit control; the invisible
   return value is unchanged.

## NOTE

```
Maintainer: 'Branimir K. Hackenberger <branimir@sciom.hr>'

New submission

Suggests or Enhances not in mainstream repositories:
  cmdstanr
Availability using Additional_repositories specification:
  cmdstanr   yes   https://stan-dev.r-universe.dev
```

**cmdstanr** is listed under `Suggests` because it is not available on
CRAN.  It is hosted on the Stan r-universe at
<https://stan-dev.r-universe.dev>.  The `Additional_repositories` field
in DESCRIPTION points to that repository.

All functions that require cmdstanr (`bdeb_fit()`) check for its
availability at runtime via `requireNamespace()` and verify that the
CmdStan toolchain is installed via `cmdstanr::cmdstan_path()`.  The
package loads, prints, and passes all 903 tests without cmdstanr
installed.

## Test coverage

* **903 unit tests** across 17 test files, covering all 26 exported
  functions.
* Test categories:
  - Data preparation and validation
  - Prior system and Stan data mapping
  - Model specification and type/data contracts
  - Fitting interface and S3 methods
  - Convergence diagnostics with bad-fit detection
  - Posterior predictive checks and predictions
  - DEBtox: EC50/NEC extraction and dose-response
  - Plotting for all 4 model types
  - Stan model file integrity
  - API contract tests
  - Scientific consistency: Arrhenius, energy conservation,
    kappa-rule, monotonicity, EC50 identity
  - Snapshot/regression: dataset stability, prior defaults,
    known numerical values
  - Deep validation: numerical stability, convergence detection,
    derived quantity formulas, dimensional analysis
* Tests use mock `bdeb_fit` objects (including bad-diagnostics mocks)
  to exercise downstream functions without requiring CmdStan.

## Test environments

* Local: Ubuntu 22.04.5 LTS, R 4.5.2, CmdStan 2.36.0
* `R CMD build` + `R CMD check --as-cran`: 0 errors, 0 warnings,
  1 note (cmdstanr in Additional_repositories)

## Downstream dependencies

None (new package).
