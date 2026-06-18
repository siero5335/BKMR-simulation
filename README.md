# Pure BKMR simulation code

This directory contains a fresh, standalone BKMR simulation code base.

The code uses only:

- `bkmr`
- `posterior`
- base R / recommended R packages

It does not source benchmark helper files from another repository and does not
load any acceleration package. Data are generated with `bkmr::SimData()`, and
model fitting is performed with `bkmr::kmbayes()`.

For binomial BKMR fits, starting values are supplied explicitly through the
official `bkmr::kmbayes(starting.values = ...)` argument:

- `h.hat = 1`
- `beta = 0`
- `sigsq.eps = 1`
- `r = 1`
- `lambda = 10`
- `delta = 1`
- `ystar = 0.5` for cases and `-0.5` for controls

This avoids relying on the internal probit GLM initialization, which can become
unstable under small-sample separation while leaving the BKMR model, priors, and
Metropolis-Hastings proposal settings otherwise unchanged.

Convergence diagnostics are computed with the `posterior` package following
Vehtari et al. (2021): `posterior::rhat()` reports the maximum of
rank-normalized split R-hat and rank-normalized folded-split R-hat, and
`posterior::ess_bulk()` and `posterior::ess_tail()` are reported separately.
The primary diagnostic flags use R-hat > 1.01 and bulk-ESS or tail-ESS < 400.

## Completed run

The 8-worker run completed at 2026-06-16 23:56:10 JST.

The completed outputs are stored under:

- `results/pure_bkmr_20260616_startvalues/`

Key files:

- `results/pure_bkmr_20260616_startvalues/all_task_results.csv`: 431 tasks, 431 ok, 0 errors
- `results/pure_bkmr_20260616_startvalues/overall_runtime.csv`: elapsed 16029.56 seconds, 8 workers
- `results/pure_bkmr_20260616_startvalues/parameter_diagnostics_by_task/`: 431 task-specific parameter
  diagnostic CSV files

The `converged_hard` and `converged_strict` columns are diagnostic flags, not
fit-completion indicators. Under the strict Vehtari-style criteria, many
completed BKMR fits are flagged for R-hat and/or ESS. These flags should be
reported as convergence-diagnostic results.

Document revision prompts for external review are stored in:

- `document_revision_prompts/`

## Scripts

- `code/run_balanced_main.R`
- `code/run_balanced_ablation.R`
- `code/run_balanced_longrun_seed3.R`
- `code/run_balanced_n40_longrun_single.R`
- `code/run_default_support.R`
- `code/run_public_verification.R`
- `code/run_all_parallel.R`
- `code/rerun_duplicate_parameter_diagnostics.R`

Set `BKMR_WORKERS=8` or another value before running `run_all_parallel.R` to
control scenario-level parallelism.

`code/rerun_duplicate_parameter_diagnostics.R` is included because some scenario
families intentionally repeat the same data set and iteration count while
changing only the chain seed set. It creates the canonical
`results/pure_bkmr_20260616_startvalues/parameter_diagnostics_by_task/`
directory with one parameter diagnostic file per completed task.
