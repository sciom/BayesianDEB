# Response to the Second-Round Review

We thank the editor and reviewer for the careful second reading of
*BayesianDEB: A Bayesian Framework for Dynamic Energy Budget Modelling
in R*.  Below we reproduce each point of the review of **22 May 2026**
(in *italics*), followed by our response and the specific files in which
the change was made.  All work refers to package version **0.2.0**; the
replication archive is `BayesianDEB_replication.zip` and the revised
manuscript is `bayesiandeb.tex`.

A one-line summary of the changes:

* a single master replication script (`replicate_all.R`) reproduces
  every manuscript number in about one minute from cached draws;
* every code listing shown in the manuscript now executes in the
  replication scripts and prints its output;
* all manuscript outputs are produced with a fixed seed from the shipped
  cache, so they reproduce bit-identically;
* the diagnostics `print` and `plot` methods are now compact;
* the `bdeb_diagnose` example family moved from `\dontrun{}` to a
  `cmdstanr`-gated `\donttest{}` block;
* the dose-response figure now draws the reference lines described in
  its caption;
* the README installation instructions no longer pull the stale CRAN
  version.

---

## 1. Package version must match the CRAN version

> *The submitted package version is not in line with the package
> published on CRAN. We require that both versions match.*

CRAN currently holds 0.1.4, whereas the manuscript and this submission
describe 0.2.0.  We have submitted 0.2.0 to CRAN so that the published
package matches the manuscript; the accompanying `cran-comments.md`
documents the version bump.  `DESCRIPTION` declares `Version: 0.2.0`.
The package passes `R CMD check --as-cran` with 0 errors / 0 warnings /
1 note (the note is the informational `cmdstanr`-in-
`Additional_repositories` line discussed below).

---

## 2. README should include cmdstanr installation instructions

> *The README file should also include instructions for the
> installation of cmdstanr.*

Added.  `README.md` (Installation section) now gives the two-step
procedure explicitly: install `cmdstanr` from the Stan r-universe and
then build the CmdStan toolchain once,

```r
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

and states that `bdeb_fit()` checks for the toolchain at runtime.

---

## 3. A single replication script that runs within one hour

> *Journal of Statistical Software requires that a single replication
> script that can be run within one hour is provided. ... parts of the
> results require that the scripts in "sbc" are run and these scripts
> are supposed to take more than 20 hours to finish. We would thus ask
> that the authors arrange their replication scripts into a smaller
> number of R files and that they also work on reducing the
> computational time needed to replicate results.*

We have added a single master script, `replication/replicate_all.R`,
that reproduces **every figure, table and printed number** in the
manuscript with one command:

```sh
Rscript replicate_all.R
```

By default it does **not** run any MCMC: it loads the archived posterior
draws shipped in `outputs/` and `sbc/` and regenerates all outputs in
**about one minute** on a laptop — well inside the one-hour budget.  The
long simulation-based-calibration runs are never executed by this
script; their archived rank matrices (`sbc/*.rds`) reproduce Figures 3-4
and Tables 4-5 directly, and regenerating them from scratch (45 min to
12 h per model) is documented separately in `sbc/README.md` and is not
part of the replication budget.  For readers who wish to recompute, a
full refit is available via `BDEB_RECOMPUTE=true BDEB_MODE=full`
(~3 h) and a quick refit via `BDEB_MODE=lite` (~30 min).

**Reference.** `replication/replicate_all.R`, `replication/README.md`,
`replication/sbc/README.md`.

---

## 4. `example("bdeb_diagnose")` should be executable (`\donttest`)

> *example("bdeb_diagnose") shows that this example is in \dontrun.
> This should, in principle, be executable code and we expect it to be
> in donttest instead so as to allow its execution using 'example',
> similarly to what has been implemented for bdeb_fit.*

Done.  `bdeb_diagnose()` and the `bdeb_diagnostics` `print`/`summary`/
`plot` methods — together with `bdeb_loo()`, `bdeb_derived()` and the
deprecated `bdeb_summary()` — now wrap their examples in a
`cmdstanr`-gated `\donttest{}` block of the same form used for
`bdeb_fit()`:

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
present.  No exported function now uses `\dontrun{}`.

**Reference.** `R/diagnostics.R`, `R/fit.R`, `R/ppc.R`, `R/debtox.R`.

---

## 5. The `plot` method for `bdeb_diagnostics` could be more readable

> *We believe that the plot method for 'bdeb_diagnostics' could be
> improved to allow a shorter and more readable display on screen.*

`plot.bdeb_diagnostics()` now hides the per-time-point latent states
(`x_sol[i,j]`, `L_hat[i]`, ...) by default and plots only the scalar
model parameters, so the R-hat / ESS panels are short and legible
instead of carrying dozens of latent rows.  A new `full = TRUE`
argument restores the complete plot, and the subtitle reports how many
latent rows were hidden.  This mirrors the behaviour already adopted for
`print.bdeb_diagnostics()`.

**Reference.** `R/diagnostics.R` (`plot.bdeb_diagnostics`); `NEWS.md`.

---

## 6. Replication output differs from the article (LOO comparison)

> *Results obtained with the replication material ("full" version) are
> still different from the results shown in the article. For instance,
> on page 13, the loo_compare for the lognormal model differs.*

The LOO comparison on this individual is genuinely Monte-Carlo
sensitive: with only `n = 13` observations a few Pareto-`k̂` values
exceed 0.7, so the importance-sampling ELPD difference fluctuates by a
fraction of its standard error between runs (we observed values in the
range -0.7 to -1.1, all with SE ~1.2).  We have:

* fixed a global seed (`seed = 42`) in every fit so the shipped result
  is deterministic;
* made `replicate_all.R` reproduce the manuscript numbers from the
  cached draws, so the printed `Δ`ELPD matches the article exactly;
* reframed the conclusion in the text as qualitative: because
  `|Δ`ELPD`|` is smaller than its standard error, neither observation
  model is decisively preferred, and the exact difference should be read
  as indicative rather than definitive.

The manuscript now prints `lognormal -1.0` (SE 1.3) and states plainly
that this difference is within its standard error.

**Reference.** `replication/00_setup.R` (seed), `replication/01_illustrations.R`,
manuscript Section 5.1.

---

## 7. README installation pulls the old CRAN version (0.1.4)

> *The instructions in the README.md suggest installing with a repos
> option that includes CRAN, which installs version 0.1.4 and not the
> submitted package version 0.2.0.*

The installation instructions have been separated.  `BayesianDEB` is now
installed with a plain `install.packages("BayesianDEB")` (which will
resolve to 0.2.0 once it is on CRAN), and `cmdstanr` is installed from
the Stan r-universe in a **separate** call (see point 2).  The combined
`repos = c("https://stan-dev.r-universe.dev", getOption("repos"))`
command that could shadow the CRAN version of `BayesianDEB` has been
removed from the README.

**Reference.** `README.md` (Installation section).

---

## 8. Replication material is incomplete (missing manuscript code)

> *Some code included in the manuscript which should be executable code
> is not contained ... E.g., on page 7, prior_species("Daphnia_magna",
> type = "debtox").*

Every code listing shown in the manuscript is now executed verbatim in
the replication scripts.  In particular,
`prior_species("Daphnia_magna", type = "debtox")` is called and printed
in `replication/01_illustrations.R`, alongside the
`prior_species("Eisenia_fetida")` listing.  We re-audited the manuscript
listing by listing and confirmed each one appears in a replication
script.

**Reference.** `replication/01_illustrations.R`.

---

## 9. `bdeb_diagnose(fit)` output differs and is too long

> *Some results obtained differ. E.g., compared to page 12,
> bdeb_diagnose(fit) reports a different number of divergent transitions
> and prints a very long table including x_sol[i,j] and L_hat[i] rows.*

Two changes address this:

* **Reproducibility.** With the fixed seed and cached draws (point 6),
  the page-12 output is reproduced exactly by `replicate_all.R`; the
  manuscript now reports a single divergent transition (1 of 4000 draws,
  < 0.1%) and discusses it.
* **Length.** `print.bdeb_diagnostics()` now shows only the scalar model
  parameters by default and hides the per-time-point latent states; the
  manuscript output on p. 12 is correspondingly compact (eight parameter
  rows, no `x_sol`/`L_hat`).  The full table remains available via
  `print(x, full = TRUE)` or `summary(x)$table`.

**Reference.** `R/diagnostics.R`, `replication/01_illustrations.R`,
manuscript Section 5.1.

---

## 10. `bdeb_derived()` result is created but never printed

> *The replication material creates der <- bdeb_derived(fit, ...) but
> does not print it, so it is not easy to compare with the manuscript.*

Fixed.  `replication/01_illustrations.R` now prints the derived-quantity
table immediately after constructing it, in the same layout as the
manuscript:

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

These warnings are now explicitly anticipated and explained.  In
`replication/01_illustrations.R` a comment precedes the `loo_compare()`
call noting that, with `n = 13`, a few `k̂ > 0.7` are expected and that
the importance-sampling ELPD is therefore indicative.  The manuscript
(Section 5.1) makes the same point in the running text: the Pareto-`k̂`
warning is reported, interpreted (the approximation is unreliable for a
few points), and used to justify reading the comparison qualitatively.

**Reference.** `replication/01_illustrations.R`, manuscript Section 5.1.

---

## 12. Figure caption refers to lines that are not in the figure

> *The caption of Figure 2 discusses something dashed red and dotted
> green. However, it is unclear to what in the figure one refers.*

`plot_dose_response()` now actually draws these reference lines: a
dashed red vertical line at the posterior-median EC$_{50}$ and a dotted
green vertical line at the posterior-median NEC, with a matching legend
subtitle.  Non-finite draws are dropped and the view is clipped with
`coord_cartesian()` so degenerate draws no longer create vertical
artefacts.  The Figure 2 caption now matches the rendered figure
element by element.

**Reference.** `R/debtox.R` (`plot_dose_response`), manuscript
Figure 2 caption.

---

## 13. Unclear which part of the replication reproduces Figures 3 and 4

> *It is unclear what part of the replication material reproduces
> Figures 3 and 4.*

`replication/README.md` now contains an explicit *Figure -> script*
table mapping every manuscript figure to the script and cached object
that produces it.  Figures 3 and 4 are listed as:

| Figure | File | Produced by |
| --- | --- | --- |
| Fig. 3 | `fig_sbc_ranks.pdf` | `02_validation.R` from `sbc/sbc_results.rds` |
| Fig. 4 | `fig_sbc_hierarchical.pdf` | `02_validation.R` from `sbc/sbc_hierarchical_results.rds` |

The same table maps all remaining figures and the corresponding
manuscript sections.

**Reference.** `replication/README.md` (Figure -> script and
Section -> script tables).

---

We thank the editor and reviewer again for the constructive feedback,
which has further improved the package and its replication material.

Branimir K. Hackenberger
on behalf of the authors
