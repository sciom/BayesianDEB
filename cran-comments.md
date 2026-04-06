## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.

### NOTEs

1. **New submission / cmdstanr not in mainstream repositories**

   ```
   Maintainer: 'Branimir K. Hackenberger <branimir@sciom.hr>'
   New submission
   Suggests or Enhances not in mainstream repositories: cmdstanr
   Availability using Additional_repositories specification:
     cmdstanr   yes   https://stan-dev.r-universe.dev
   ```

   **cmdstanr** is listed under `Suggests` because it is not available
   on CRAN.  It is hosted on the Stan r-universe at
   <https://stan-dev.r-universe.dev>.  The `Additional_repositories`
   field in DESCRIPTION points to that repository.

   All functions that require cmdstanr (`bdeb_fit()`) check for its
   availability at runtime via `requireNamespace()` and verify that
   the CmdStan toolchain is installed via `cmdstanr::cmdstan_path()`.
   The package loads, prints, and passes **all 753 tests** without
   cmdstanr installed.

2. **unable to verify current time** — transient network/clock issue
   on the test machine, not related to the package.

## Test coverage

* **753 unit tests** across 14 test files, covering all 25 exported
  functions.
* Test categories:
  - Data preparation and validation (74 tests)
  - Prior system and Stan data mapping (145 tests)
  - Model specification and type/data contracts (73 tests)
  - Fitting interface and S3 methods (14 tests)
  - Convergence diagnostics with bad-fit detection (47 tests)
  - Posterior predictive checks and predictions (29 tests)
  - DEBtox: EC50/NEC extraction and dose-response (25 tests)
  - Plotting for all 4 model types (17 tests)
  - Stan model file integrity (35 tests)
  - API contract tests (36 tests)
  - Scientific consistency: Arrhenius, energy conservation,
    kappa-rule, monotonicity, EC50 identity (41 tests)
  - Snapshot/regression: dataset stability, prior defaults,
    known numerical values (64 tests)
  - Deep validation: numerical stability, convergence detection,
    derived quantity formulas, dimensional analysis (101 tests)
* Tests use mock `bdeb_fit` objects (including bad-diagnostics mocks)
  to exercise downstream functions without requiring CmdStan.

## Test environments

* Local: Ubuntu 22.04.5 LTS, R 4.5.2, CmdStan 2.36.0
* `R CMD build` + `R CMD check --as-cran`: 0 errors, 0 warnings,
  2 notes (new submission + cmdstanr)

## Downstream dependencies

None (new package).
