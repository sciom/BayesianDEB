# cran-comments

## Submission comment (as sent to CRAN)

This is a new version (0.2.1) superseding the current CRAN release 0.1.4.
It accompanies a Journal of Statistical Software submission, which requires
the published package version to match the manuscript.

### R CMD check results
0 errors | 0 warnings | 1 note
Tested on: Ubuntu 22.04.5 LTS, R 4.5.2 (local R CMD check --as-cran).

### The one NOTE
```
Suggests or Enhances not in mainstream repositories: cmdstanr
Availability using Additional_repositories specification:
  cmdstanr   yes   https://stan-dev.r-universe.dev
```

cmdstanr is in Suggests (not Imports) because it is not on CRAN; it is
hosted on the Stan r-universe, to which the Additional_repositories field
in DESCRIPTION points. The package installs, loads, and passes its full
test suite without cmdstanr. Only model fitting (`bdeb_fit()`) needs it;
that function checks availability at runtime via `requireNamespace()` and
verifies the CmdStan toolchain via `cmdstanr::cmdstan_path()`, returning an
informative message if either is missing. All examples using cmdstanr are
wrapped in toolchain-gated `\donttest{}` blocks.

### Downstream dependencies
None. The package has no reverse dependencies on CRAN; the existing
0.1.4 reverse-dependency set is empty, so this update affects nothing.

---

## What changed in 0.2.1

Version 0.2.0 added a full S3 class system (print/summary/plot for every
exported class), additional `lm`-style methods on `bdeb_fit` (`confint`,
`fitted`, `residuals`, `nobs`, `vcov`, `logLik`), LOO/WAIC comparison, and
an expanded test suite.  Version 0.2.1 (this submission) adds the
second-round JSS revision: all exported examples are now cmdstanr-gated
`\donttest{}` blocks (no `\dontrun{}` remains), the `bdeb_diagnostics`
print and plot methods hide per-time-point latent states by default for a
compact display, and the replication material ships a single
`replicate_all.R` reproducing every manuscript result in about one minute.
See `NEWS.md` for details.

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
   As of 0.2.1 the remaining `bdeb_fit()`-style examples (which require
   the external CmdStan toolchain) are wrapped in cmdstanr-gated
   `\donttest{}` blocks rather than `\dontrun{}`.

3. **Suppressible output.**  Direct `print()` calls inside
   `bdeb_diagnose()` and `bdeb_ec50()` were removed; all user-facing
   output now goes through `cli` / `message()` and can be silenced
   with `suppressMessages()`.  Both functions also gained a new
   `verbose = TRUE` argument for explicit control; the invisible
   return value is unchanged.

## Test coverage

* **903 unit tests** across 17 test files, covering all 26 exported
  functions.
* Test categories: data preparation and validation; prior system and
  Stan data mapping; model specification and type/data contracts;
  fitting interface and S3 methods; convergence diagnostics with
  bad-fit detection; posterior predictive checks and predictions;
  DEBtox EC50/NEC extraction and dose-response; plotting for all 4
  model types; Stan model file integrity; API contract tests;
  scientific consistency (Arrhenius, energy conservation, kappa-rule,
  monotonicity, EC50 identity); snapshot/regression; deep validation
  (numerical stability, convergence detection, derived-quantity
  formulas, dimensional analysis).
* Tests use mock `bdeb_fit` objects (including bad-diagnostics mocks)
  to exercise downstream functions without requiring CmdStan.

## Test environments

* Local: Ubuntu 22.04.5 LTS, R 4.5.2, CmdStan 2.36.0
* `R CMD build` + `R CMD check --as-cran`: 0 errors, 0 warnings,
  1 note (cmdstanr in Additional_repositories)
