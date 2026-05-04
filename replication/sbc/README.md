# SBC scripts (Section 6.1)

These scripts perform Simulation-Based Calibration (Talts et al. 2020) for
the four BayesianDEB model types.  They are **long-running** and are kept
out of the main replication flow (`01_`-`03_`); the main flow loads the
pre-computed RDS results in this directory instead.

| Script                  | Model        | Replications | Approx. runtime (16-core) |
|-------------------------|--------------|--------------|---------------------------|
| `sbc_individual.R`      | individual   | 500          | ~45 min                   |
| `sbc_hierarchical.R`    | hierarchical | 50           | ~5 h                      |
| `sbc_growth_repro.R`    | growth_repro | 100          | ~2-4 h                    |
| `sbc_debtox.R`          | debtox       | 100          | ~6-12 h                   |

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
