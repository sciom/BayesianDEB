# SBC scripts (Section 6.1)

These scripts perform Simulation-Based Calibration (Talts et al. 2020) for
the BayesianDEB model types.  They are **long-running** and are **not run**
by `replicate_all.R` or the numbered scripts.

**You do not need to run anything in this folder to reproduce the
manuscript.**  The paper reports SBC only for the *individual* model
(Figure 3 / Table 4) and the *hierarchical* model (Figure 4 / Table 5);
the archived rank matrices for both are shipped here as `.rds`.
`02_validation.R` loads them and rebuilds those figures and tables in
seconds.

| Script                  | Model        | Reps | Runtime (16-core) | In manuscript?     |
|-------------------------|--------------|------|-------------------|--------------------|
| `sbc_individual.R`      | individual   | 500  | ~45 min           | Yes (Fig 3, Tab 4) |
| `sbc_hierarchical.R`    | hierarchical | 50   | ~5 h              | Yes (Fig 4, Tab 5) |
| `sbc_growth_repro.R`    | growth_repro | 100  | ~2-4 h            | No (exploratory)   |
| `sbc_debtox.R`          | debtox       | 100  | ~6-12 h           | No (exploratory)   |

The `growth_repro` and `debtox` SBC scripts are exploratory, provided only
for completeness; their results are not shown in the paper and are not
required by any part of the replication.

## Pre-computed results

The following RDS files are bundled so the main scripts can summarise SBC
without rerunning:

* `sbc_results.rds`               — output of `sbc_individual.R`
* `sbc_hierarchical_results.rds`  — output of `sbc_hierarchical.R`

## Re-running

Each script is self-contained.  To regenerate the RDS:

```bash
Rscript replication/sbc/sbc_individual.R       # ~45 min
Rscript replication/sbc/sbc_hierarchical.R     # ~5 h
Rscript replication/sbc/sbc_growth_repro.R     # ~2-4 h
Rscript replication/sbc/sbc_debtox.R           # ~6-12 h
```

The output RDS is written next to the script (e.g.
`sbc_results.rds`).
