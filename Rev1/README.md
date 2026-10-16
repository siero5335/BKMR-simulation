# First-revision additional analyses

Prepared for public code and numerical-result sharing on 16 October 2026. This directory supplements the original `code/` and `results/` at the repository root. It does not change the original 431-task design, the six-parameter primary criterion, or the original classifications.

## Contents and manuscript correspondence

| Supplementary Information | Purpose | Published numerical results | Rerun stages |
|---|---|---|---|
| Section 3 | Postprocessed h at all observed exposure vectors, with est.h=FALSE and varsel=FALSE | `h_observed_points/` | `beta`, then `h`, then `h-summary` |
| Section 4 | Beta bias, RMSE, interval coverage, and MCSE | `tables/beta_performance.csv`, `tables/beta_mcse.csv`, `coauthor_additions/beta_estimand_by_task.csv` | `audit`, `beta`, `beta-summary` |
| Section 5, Figure 5.1, Table 5.2 | Sensitivity to initial states and iteration budgets | `initialization_extension/`, `coauthor_additions/` | `beta`, `init-inputs`, `init`, `init-summary`, `figures` |
| Section 7 | Numerical non-completion and matrix checks | `failure_case/` | `audit`, `beta --task 174`, `failure-trace` |
| Section 10 | Parameter-flag overlap and lambda acceptance summaries | `tables/parameter_overlap.csv`, `tables/lambda_acceptance.csv` | `audit` |
| Section 11 | Reporting template and decision workflow | Provided in the SI; no additional fitting experiment | None |

`coauthor_additions/` retains its historical folder name to preserve file paths. It contains numerical checks added during coauthor review, including the descriptive endpoint-MCSE counts and the selected long-chain reference comparisons. These are dependent comparisons, not new independent simulation replicates.

## Quick verification: no MCMC

From the repository root:

```sh
python3 Rev1/code/verify_results.py
```

Requirements: Python 3.9+ (standard library only) and base R on PATH. If necessary, set `RSCRIPT` to the Rscript executable path. The command uses temporary output directories and checks:

- the coauthor summary JSON and both supporting comparison tables;
- four initialization/budget/warmup comparison tables;
- all 29,680 point-level h diagnostic records and their task-level primary classifications;
- all 430 completed-task primary classifications against the original summary records.

These checks reaggregate stored numerical results. They do **not** independently recompute R-hat or ESS from historical draws.

## Execution code

- `code/manage.py`: public entry point; prepares a new isolated workspace and runs explicit stages.
- `code/runtime/`: public rerun adapters that preserve the original statistical calculations while replacing dependencies on unavailable local checkpoints and execution paths.
- `code/original/`: unchanged historical analytical scripts, retained for method inspection and provenance. These are not the recommended public entry points and can refer to omitted historical checkpoints.
- `code/SOURCE_INDEX.csv`: source-relative paths and SHA-256 checksums for the historical scripts.
- `reference/`: consolidated original parameter diagnostics, core diagnostics associated with the h analysis, and twelve exact simulated-data/start specifications for the initialization experiment.
- `environment/`: recorded session information and required package versions.

See [REPRODUCING.md](REPRODUCING.md) for commands, dependencies, expected outputs, and computational stages. See [PROVENANCE.md](PROVENANCE.md) for the distinction between historical verification and a new rerun.

## Scope and interpretation

The beta estimand is the covariate's probit coefficient, with generating value 0.5. It is not the overall mixture effect.

The h analysis covers 200 planned datasets (n=40, 80, 160, 320), of which 199 completed and one had the original numerical non-completion. h was sampled after fitting with est.h=FALSE, with chain and retained-iteration alignment preserved. These diagnostics concern h at observed exposure vectors; they do not establish convergence of every posterior quantity. `h_summary.csv` contains both all-point and nested-five summaries; do not pool those rows.

The initialization extension uses 12 deliberately selected datasets, three initial-state sets, four chains per set, and prefixes at 10k/20k/40k/100k iterations. The same sampler seeds are used across initialization sets and the prefixes reuse trajectories. These comparisons are descriptive and dependent. Diagnostic non-pass alone does not establish weak identification, and small MCSE does not establish adequate exploration.

## Included and omitted data

Included: the original lightweight summary CSVs, all-point h diagnostic CSV (29,680 rows), seed/registry/plan information, exact simulated inputs and initial states for 12 extension cases, execution code, and environment records. Existing numerical results are unchanged.

Omitted: raw draws, fitted models, chain RDS files, checkpoints, intermediate caches, and manuscript/response drafts. The input specification is text (`reference/initialization_specs.R`), exported using hexadecimal floating-point values to preserve the exact inputs and starts; it contains no MCMC trajectories. Full reruns regenerate needed data and baseline draws locally. The runtime does not require access to the author's original computer or private files.

`.gitignore` excludes local rerun outputs. Do not force-add generated RDS/checkpoint files.
