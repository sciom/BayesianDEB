# Response to the Second-Round Review

We are grateful to the editor and reviewer for reading *BayesianDEB: A
Bayesian Framework for Dynamic Energy Budget Modelling in R* so closely a
second time.  Each comment from the review of **22 May 2026** is quoted
below (in *italics*), with our reply underneath and the files we touched
to address it.  Everything here refers to package version **0.2.1**; the
replication archive is `BayesianDEB_replication.zip` and the revised
manuscript is `bayesiandeb.tex`.

In short:

* a single master replication script (`replicate_all.R`) rebuilds every
  manuscript number in about a minute from cached draws;
* every code listing in the manuscript now runs, and prints its output,
  in the replication scripts;
* all manuscript outputs come from the shipped cache under a fixed seed,
  so they reproduce bit for bit;
* the diagnostics `print` and `plot` methods are now compact;
* the `bdeb_diagnose` examples have moved from `\dontrun{}` to a
  `cmdstanr`-gated `\donttest{}` block;
* the dose-response figure now draws the reference lines its caption
  describes;
* the README installation instructions no longer pull the stale CRAN
  version.

---

## 1. Package version must match the CRAN version

> *The submitted package version is not in line with the package
> published on CRAN. We require that both versions match.*

When the reviewer wrote, CRAN still held 0.1.4 while the manuscript and
this submission describe 0.2.1.  That gap is now closed: 0.2.1 went up on
CRAN on 2026-06-17, so the published package and the manuscript describe
the same release.  The version bump is documented in `cran-comments.md`,
and `DESCRIPTION` declares `Version: 0.2.1`.  `R CMD check --as-cran`
reports 0 errors, 0 warnings and a single note (the informational
`cmdstanr`-in-`Additional_repositories` line, discussed below).

---

## 2. README should include cmdstanr installation instructions

> *The README file should also include instructions for the
> installation of cmdstanr.*

Added.  The Installation section of `README.md` now spells out the two
steps: first install `cmdstanr` from the Stan r-universe, then build the
CmdStan toolchain once,

```r
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

and the text notes that `bdeb_fit()` checks for the toolchain at run
time.

---

## 3. A single replication script that runs within one hour

> *Journal of Statistical Software requires that a single replication
> script that can be run within one hour is provided. ... parts of the
> results require that the scripts in "sbc" are run and these scripts
> are supposed to take more than 20 hours to finish. We would thus ask
> that the authors arrange their replication scripts into a smaller
> number of R files and that they also work on reducing the
> computational time needed to replicate results.*

There is now a single master script, `replication/replicate_all.R`, that
rebuilds every figure, table and printed number in the manuscript from
one command:

```sh
Rscript replicate_all.R
```

Out of the box it runs no MCMC at all: it reads the archived posterior
draws in `outputs/` and `sbc/` and regenerates everything in roughly a
minute on a laptop, comfortably under the one-hour limit.  The lengthy
simulation-based-calibration runs stay out of this path.  Their rank
matrices are shipped (`sbc/*.rds`) and reproduce Figures 3-4 and
Tables 4-5 on their own; recomputing them from scratch takes anywhere
from 45 min to 12 h per model and is described separately in
`sbc/README.md`, outside the replication budget.  Anyone who does want to
refit can set `BDEB_RECOMPUTE=true BDEB_MODE=full` for a complete run
(~3 h) or `BDEB_MODE=lite` for a quick one (~30 min).

**Reference.** `replication/replicate_all.R`, `replication/README.md`,
`replication/sbc/README.md`.

---

## 4. `example("bdeb_diagnose")` should be executable (`\donttest`)

> *example("bdeb_diagnose") shows that this example is in \dontrun.
> This should, in principle, be executable code and we expect it to be
> in donttest instead so as to allow its execution using 'example',
> similarly to what has been implemented for bdeb_fit.*

Done.  The examples for `bdeb_diagnose()`, the `bdeb_diagnostics`
`print`/`summary`/`plot` methods, `bdeb_loo()`, `bdeb_derived()` and the
deprecated `bdeb_summary()` are now wrapped in the same `cmdstanr`-gated
`\donttest{}` block we already use for `bdeb_fit()`:

```r
\donttest{
if (requireNamespace("cmdstanr", quietly = TRUE) &&
    nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
  data(eisenia_growth)
  dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
  fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
                  chains = 2, iter_warmup = 200, iter_sampling = 200,
                  refresh = 0)
  print(bdeb_diagnose(fit))
}
}
```

so `example("bdeb_diagnose")` runs whenever a CmdStan toolchain is
available.  None of the exported functions rely on `\dontrun{}` any more.

**Reference.** `R/diagnostics.R`, `R/fit.R`, `R/ppc.R`, `R/debtox.R`.

---

## 5. The `plot` method for `bdeb_diagnostics` could be more readable

> *We believe that the plot method for 'bdeb_diagnostics' could be
> improved to allow a shorter and more readable display on screen.*

`plot.bdeb_diagnostics()` no longer crowds the display with
per-time-point latent states (`x_sol[i,j]`, `L_hat[i]`, and so on).  By
default it shows only the scalar model parameters, so the R-hat and ESS
panels stay short and easy to read rather than running to dozens of
latent rows.  Passing `full = TRUE` brings back the complete plot, and
the subtitle notes how many latent rows were left out.  This follows what
`print.bdeb_diagnostics()` already does.

**Reference.** `R/diagnostics.R` (`plot.bdeb_diagnostics`); `NEWS.md`.

---

## 6. Replication output differs from the article (LOO comparison)

> *Results obtained with the replication material ("full" version) are
> still different from the results shown in the article. For instance,
> on page 13, the loo_compare for the lognormal model differs.*

This comparison really is sensitive to Monte-Carlo noise on this
individual.  With only `n = 13` observations, a handful of Pareto-`k̂`
values climb above 0.7, so the importance-sampling ELPD difference moves
by a fraction of its standard error from run to run; across our reruns it
landed around -1 every time, always with SE ~1.3 and always within one
standard error of zero.  We made three changes:

* fixed a global seed (`seed = 42`) in every fit, so the shipped result
  is deterministic;
* had `replicate_all.R` reproduce the manuscript numbers from the cached
  draws, so the printed `Δ`ELPD matches the article exactly;
* rewrote the conclusion as a qualitative one: since `|Δ`ELPD`|` is
  smaller than its standard error, neither observation model wins
  clearly, and the precise difference is best read as indicative.

The manuscript prints `lognormal -1.0` (SE 1.3) and says outright that
the difference sits within its standard error.

**Reference.** `replication/00_setup.R` (seed), `replication/01_illustrations.R`,
manuscript Section 5.1.

---

## 7. README installation pulls the old CRAN version (0.1.4)

> *The instructions in the README.md suggest installing with a repos
> option that includes CRAN, which installs version 0.1.4 and not the
> submitted package version 0.2.0.*

We have split the installation instructions.  `BayesianDEB` is now
installed with a plain `install.packages("BayesianDEB")`, which now
resolves to 0.2.1 from CRAN, while `cmdstanr` is installed separately
from the Stan r-universe (see point 2).  The earlier combined
`repos = c("https://stan-dev.r-universe.dev", getOption("repos"))` call,
which could have masked the CRAN copy of `BayesianDEB`, is gone from the
README.

**Reference.** `README.md` (Installation section).

---

## 8. Replication material is incomplete (missing manuscript code)

> *Some code included in the manuscript which should be executable code
> is not contained ... E.g., on page 7, prior_species("Daphnia_magna",
> type = "debtox").*

Every code listing in the manuscript is now run in the replication
scripts.  `prior_species("Daphnia_magna", type = "debtox")`, for example,
is called and printed in `replication/01_illustrations.R` next to the
`prior_species("Eisenia_fetida")` listing.  We went back through the
manuscript one listing at a time and checked that each appears in a
replication script.

**Reference.** `replication/01_illustrations.R`.

---

## 9. `bdeb_diagnose(fit)` output differs and is too long

> *Some results obtained differ. E.g., compared to page 12,
> bdeb_diagnose(fit) reports a different number of divergent transitions
> and prints a very long table including x_sol[i,j] and L_hat[i] rows.*

Two things address this:

* **Reproducibility.** With the fixed seed and cached draws (point 6),
  `replicate_all.R` reproduces the page-12 output exactly; the manuscript
  now reports a single divergent transition (1 of 4000 draws, < 0.1%) and
  comments on it.
* **Length.** `print.bdeb_diagnostics()` now defaults to the scalar model
  parameters and hides the per-time-point latent states, so the p. 12
  output is compact (eight parameter rows, no `x_sol`/`L_hat`).  The full
  table remains available through `print(x, full = TRUE)` or
  `summary(x)$table`.

**Reference.** `R/diagnostics.R`, `replication/01_illustrations.R`,
manuscript Section 5.1.

---

## 10. `bdeb_derived()` result is created but never printed

> *The replication material creates der <- bdeb_derived(fit, ...) but
> does not print it, so it is not easy to compare with the manuscript.*

Fixed.  `replication/01_illustrations.R` now prints the derived-quantity
table right after it is built, laid out as in the manuscript:

```r
der <- bdeb_derived(fit,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"), f = 1.0)
print(as.data.frame(summarise_draws(der, "mean", "sd",
  "q5" = ~quantile(.x, 0.05), "q95" = ~quantile(.x, 0.95))), digits = 3)
```

**Reference.** `replication/01_illustrations.R`.

---

## 11. LOO triggers Pareto-`k̂` warnings that are not commented on

> *Warnings are triggered (Some Pareto k diagnostic values are too
> high) where one might also point them out and comment on them.*

These warnings are now flagged and explained rather than left to surprise
the reader.  In `replication/01_illustrations.R` a comment ahead of the
`loo_compare()` call points out that, at `n = 13`, a few `k̂ > 0.7` are
to be expected and the importance-sampling ELPD is therefore only
indicative.  The manuscript makes the same point in the text of
Section 5.1: it reports the Pareto-`k̂` warning, explains what it means
(the approximation is unreliable at a few points), and uses it to justify
reading the comparison qualitatively.

**Reference.** `replication/01_illustrations.R`, manuscript Section 5.1.

---

## 12. Figure caption refers to lines that are not in the figure

> *The caption of Figure 2 discusses something dashed red and dotted
> green. However, it is unclear to what in the figure one refers.*

`plot_dose_response()` now draws the lines the caption promises: a dashed
red vertical line at the posterior-median EC$_{50}$ and a dotted green
one at the posterior-median NEC, with a matching legend subtitle.
Non-finite draws are discarded and the view is clipped with
`coord_cartesian()`, so degenerate draws no longer leave vertical
artefacts.  The Figure 2 caption and the rendered figure now agree.

**Reference.** `R/debtox.R` (`plot_dose_response`), manuscript
Figure 2 caption.

---

## 13. Unclear which part of the replication reproduces Figures 3 and 4

> *It is unclear what part of the replication material reproduces
> Figures 3 and 4.*

`replication/README.md` now carries a *Figure -> script* table that ties
each manuscript figure to the script and cached object behind it.  For
the two in question:

| Figure | File | Produced by |
| --- | --- | --- |
| Fig. 3 | `fig_sbc_ranks.pdf` | `02_validation.R` from `sbc/sbc_results.rds` |
| Fig. 4 | `fig_sbc_hierarchical.pdf` | `02_validation.R` from `sbc/sbc_hierarchical_results.rds` |

The same table covers the remaining figures and the matching manuscript
sections.

**Reference.** `replication/README.md` (Figure -> script and
Section -> script tables).

---

We thank the editor and reviewer once more for feedback that has made the
package and its replication material better.

Branimir K. Hackenberger
on behalf of the authors
