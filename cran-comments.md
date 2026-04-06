## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes on package dependencies

* **cmdstanr** is listed under `Suggests` because it is not available on
  CRAN (it is hosted on the Stan r-universe at
  <https://stan-dev.r-universe.dev>). The `Additional_repositories` field
  in DESCRIPTION points to that repository so that
  `install.packages()` can find it.

* All functions that require cmdstanr (`bdeb_fit()`) check for its
  availability at runtime via `requireNamespace()` and provide an
  informative error message with installation instructions if it is
  missing. The package loads, prints, and passes all non-fitting tests
  without cmdstanr installed.

* **CmdStan** (the C++ Stan compiler) is listed under
  `SystemRequirements`. It is only needed to compile and run Stan models
  via `bdeb_fit()`.

## Test environments

* Local: Ubuntu 22.04 LTS, R 4.5.2
* R CMD check --as-cran: 0 errors, 0 warnings, 1 note (new submission)

## Downstream dependencies

None (new package).
