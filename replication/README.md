# BayesianDEB — JSS replication material

This archive reproduces all numerical results, tables and figures in

> Hackenberger BK, Djerdj T, Hackenberger DK (2026).
> *BayesianDEB: A Bayesian Framework for Dynamic Energy Budget
> Modelling in R.*  Submitted to the Journal of Statistical Software.

## Quick start

```bash
# install the package (Suggests cmdstanr)
R -q -e 'install.packages("BayesianDEB",
                          repos = c("https://stan-dev.r-universe.dev",
                                    getOption("repos")))'
R -q -e 'cmdstanr::install_cmdstan()'

# lite mode (~30 min total on a 4-core machine)
cd replication
Rscript 01_illustrations.R
Rscript 02_validation.R
Rscript 03_case_studies.R
```

To reproduce the **publication-grade** numbers in the manuscript, set

```bash
export BDEB_MODE=full        # full mode: ~3 h on 8 cores
Rscript 01_illustrations.R
Rscript 02_validation.R
Rscript 03_case_studies.R
```

| `BDEB_MODE` | Chains | Warmup | Sampling | Approx. total runtime |
|-------------|--------|--------|----------|-----------------------|
| `lite` (default) | 2 | 300 | 300 | **~30 min** (4 cores) |
| `full` | 4 | 1000 | 1000 | **~3 h** (8 cores) |

**Note:** `lite` mode is intended only as a quick API/numerics
sanity check.  R-hat / divergence diagnostics are expected to
fail in `lite` for the larger models (hierarchical §5.3, full
DEBtox §7.2); use `BDEB_MODE=full` for publication-grade results.

All output figures and tables are written to `replication/outputs/`.
The reproducibility report (BayesianDEB / cmdstanr / CmdStan / R
versions, Stan model hashes) is written to
`replication/outputs/sessionInfo.txt`.


## Mapping manuscript ↔ scripts

| Manuscript                              | Script                                 |
|-----------------------------------------|----------------------------------------|
| §5.1 Individual growth (E. fetida)      | `01_illustrations.R` (`Section 5.1`)   |
| §5.2 Growth + reproduction (F. candida) | `01_illustrations.R` (`Section 5.2`)   |
| §5.3 Hierarchical model                 | `01_illustrations.R` (`Section 5.3`)   |
| §5.4 DEBtox (eisenia_cd)                | `01_illustrations.R` (`Section 5.4`)   |
| §6.1 Simulation-based calibration       | `02_validation.R` + `sbc/` (auxiliary) |
| §6.2 Identifiability / contraction      | `02_validation.R` (Table 5, fig)       |
| §7.1 Real Eisenia data (Neuhauser 1980) | `03_case_studies.R` (`Section 7.1`)    |
| §7.2 Real DEBtox (Van Gestel 1991)      | `03_case_studies.R` (full mode only)   |
| §9   Comparison with deBInfer           | `comparison_debinfer.R`                |

The page-7 listing `prior_species("Eisenia_fetida")` is included
explicitly in `01_illustrations.R` (Section 5.1).


## Directory structure

```
replication/
├── README.md                   ← this file
├── 00_setup.R                  ← libraries, palette, BDEB_MODE
├── 01_illustrations.R          ← Section 5
├── 02_validation.R             ← Section 6
├── 03_case_studies.R           ← Section 7
├── comparison_debinfer.R       ← Section 9 (auxiliary)
├── data/
│   └── curves.txt              ← EGrowth (Mathieu, 2018) — bundled
├── sbc/
│   ├── README.md
│   ├── sbc_individual.R        ← long-running SBC runners
│   ├── sbc_hierarchical.R
│   ├── sbc_growth_repro.R
│   ├── sbc_debtox.R
│   ├── sbc_results.rds         ← pre-computed (loaded by 02_validation.R)
│   └── sbc_hierarchical_results.rds
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
