# BayesianDEB — JSS replication material

This archive reproduces all numerical results, tables and figures in

> Hackenberger BK, Djerdj T, Hackenberger DK (2026).
> *BayesianDEB: A Bayesian Framework for Dynamic Energy Budget
> Modelling in R.*  Submitted to the Journal of Statistical Software.

## Quick start

A **single script** reproduces every figure, table and printed number in
the manuscript:

```bash
# 1. Install BayesianDEB (>= 0.2.0)
R -q -e 'install.packages("BayesianDEB")'

# 2. Reproduce everything (a few minutes; no MCMC -- see below)
cd replication
Rscript replicate_all.R
```

By default `replicate_all.R` **loads the archived posterior draws** that
are shipped in `outputs/` and `sbc/`, and regenerates all figures and
tables from them.  No MCMC is run, so the whole manuscript is reproduced
**in a few minutes** (about 1 min on a laptop), well inside the JSS
one-hour budget, and the output is bit-identical to the manuscript on
every platform (Hamiltonian Monte Carlo is not bit-reproducible across
machines, so loading the archived draws is what guarantees identical
numbers).

To **refit every model from scratch** instead (requires a working
CmdStan toolchain):

```bash
R -q -e 'cmdstanr::install_cmdstan()'          # one-time
BDEB_RECOMPUTE=true BDEB_MODE=full Rscript replicate_all.R
```

| Invocation | What runs | Approx. runtime |
|------------|-----------|-----------------|
| `Rscript replicate_all.R` (default) | load archived draws, rebuild all figures/tables | **~1 min** |
| `BDEB_RECOMPUTE=true BDEB_MODE=full ...` | refit all models (4 chains, 1000+1000) | **~3 h** (8 cores) |
| `BDEB_RECOMPUTE=true BDEB_MODE=lite ...` | quick refit (2 chains, 300+300) | **~30 min** |

When refitting, `BDEB_MODE=lite` is only a fast sampler smoke-test;
R-hat / divergence diagnostics are expected to be poor for the larger
models.  Use `BDEB_MODE=full` for publication-grade refits.

The long simulation-based-calibration (SBC) runs are **never** executed
by `replicate_all.R`.  Their archived rank matrices (`sbc/*.rds`)
reproduce Figures 3-4 and Tables 4-5 in seconds; regenerating the rds
from scratch (45 min -- 12 h per model) is documented in `sbc/README.md`
and is *not* part of the replication budget.

All output figures and tables are written to `replication/outputs/`.
The reproducibility report (package / cmdstanr / CmdStan / R versions)
is written to `replication/outputs/sessionInfo.txt`.

You can still run the three numbered scripts individually
(`Rscript 01_illustrations.R`, etc.); each is cache-aware and honours
the same `BDEB_RECOMPUTE` / `BDEB_MODE` variables.


## Mapping manuscript ↔ scripts

| Manuscript                              | Script                                 |
|-----------------------------------------|----------------------------------------|
| §5.1 Individual growth (E. fetida)      | `01_illustrations.R` (`Section 5.1`)   |
| §5.2 Growth + reproduction (F. candida) | `01_illustrations.R` (`Section 5.2`)   |
| §5.3 Hierarchical model                 | `01_illustrations.R` (`Section 5.3`)   |
| §5.4 DEBtox (eisenia_cd)                | `01_illustrations.R` (`Section 5.4`)   |
| §6.1 Simulation-based calibration       | `02_validation.R` (loads `sbc/*.rds`)  |
| §6.2 Identifiability / contraction      | `02_validation.R` (Table 5, fig)       |
| §7.1 Real Eisenia data (Neuhauser 1980) | `03_case_studies.R` (`Section 7.1`)    |
| §7.2 Real DEBtox (Van Gestel 1991)      | `03_case_studies.R` (from archived fit)|
| §9   Comparison with deBInfer           | `comparison_debinfer.R`                |

### Figure → script

| Figure | File | Produced by |
|--------|------|-------------|
| Fig. 1 | `fig_debtox_rawdata.pdf`      | `01_illustrations.R` (§5.4) |
| Fig. 2 | `fig_debtox_doseresponse.pdf` | `01_illustrations.R` (§5.4, `plot_dose_response`) |
| Fig. 3 | `fig_sbc_ranks.pdf`           | `02_validation.R` from `sbc/sbc_results.rds` |
| Fig. 4 | `fig_sbc_hierarchical.pdf`    | `02_validation.R` from `sbc/sbc_hierarchical_results.rds` |
| Fig. 5 | `fig_contraction.pdf`         | `02_validation.R` (§6.2) |
| Fig. 6 | `fig_trajectory.pdf`          | `03_case_studies.R` (§7.1) |
| Fig. 7 | `fig_ppc.pdf`                 | `03_case_studies.R` (§7.1) |
| Fig. 8 | `fig_pairs.pdf`               | `03_case_studies.R` (§7.1, needs `gridExtra`) |
| Fig. 9 | `fig_trace.pdf`               | `03_case_studies.R` (§7.1) |
| Fig. 10| `fig_posterior.pdf`           | `03_case_studies.R` (§7.1) |
| Fig. 11| `fig_prior_posterior.pdf`     | `03_case_studies.R` (§7.1) |

The page-7 listings `prior_species("Eisenia_fetida")` and
`prior_species("Daphnia_magna", type = "debtox")` are both executed
explicitly in `01_illustrations.R` (Sections 5.1 and 5.4).


## Directory structure

```
replication/
├── README.md                   ← this file
├── replicate_all.R             ← SINGLE entry point (run this)
├── 00_setup.R                  ← libraries, palette, BDEB_MODE/RECOMPUTE
├── 01_illustrations.R          ← Section 5 (cache-aware)
├── 02_validation.R             ← Section 6 (SBC figs/tables from cache)
├── 03_case_studies.R           ← Section 7 (cache-aware)
├── comparison_debinfer.R       ← Section 9 (auxiliary)
├── data/
│   └── curves.txt              ← EGrowth (Mathieu, 2018) — bundled
├── sbc/
│   ├── README.md
│   ├── sbc_individual.R        ← long-running SBC runners (optional)
│   ├── sbc_hierarchical.R
│   ├── sbc_growth_repro.R      ← exploratory, NOT used in manuscript
│   ├── sbc_debtox.R            ← exploratory, NOT used in manuscript
│   ├── sbc_results.rds         ← Fig 3 / Table 4 (loaded by 02)
│   └── sbc_hierarchical_results.rds  ← Fig 4 / Table 5 (loaded by 02)
└── outputs/                    ← generated PDFs, RDS, sessionInfo.txt
```


## Bundled data

### `data/curves.txt`

Source: <https://github.com/JeromeMathieuEcology/EGrowth>
(file `curves.txt` at HEAD as of 2026-05-05; commit hash recorded
below).  Contains 16 003 records × 4 columns
(`CURVE_ID`, `time`, `bm`, `SE`).

* SHA-256: `19e131b80e8d426c5647751e056ecdb1761a6d0f85d24a54f3b34ff47cf4289b`
* Size: 326 177 bytes
* License: see EGrowth repository (CC-BY 4.0; reproduced here per
  attribution).

The file is bundled to remove the dependency on an external URL
(reviewer comment, 2026-04 round).

### Other data

The simulated datasets (`eisenia_growth`, `debtox_growth`,
`folsomia_repro`, `eisenia_cd`, `eisenia_neuhauser`) are loaded
from the installed `BayesianDEB` package via `data()`.


## Computing environment

The numbers in the manuscript were generated with the configuration
written to `outputs/sessionInfo.txt` after the last `full` run.
At submission time:

* R 4.5.x
* BayesianDEB 0.2.0
* cmdstanr 0.8.x with CmdStan 2.36.x
* OS: Ubuntu 22.04 (x86_64)
* 8-core CPU, 32 GB RAM

Slight numerical differences (≤ 1 % in posterior summaries, ≤ 4th
decimal in rates) between repeated runs reflect MCMC stochasticity
even with a fixed seed and are expected.


## Reproducibility checklist

* `set.seed(20260418)` is the global RNG seed (in `00_setup.R`).
* Each `bdeb_fit()` call uses an explicit `seed = ...`.
* Stan model source hashes are recorded in
  `outputs/sessionInfo.txt` via `bdeb_session_info()`.
* `BDEB_MODE` in `Sys.getenv()` controls runtime budget without
  editing scripts.

If a reviewer cannot reproduce a number, please open an issue at
<https://github.com/sciom/BayesianDEB/issues> with the failing
output.
