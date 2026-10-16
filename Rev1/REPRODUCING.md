# Reproduction guide

Run the commands below from the repository root. Python 3.9+ and Rscript must be available on PATH; `RSCRIPT=/path/to/Rscript` can select another executable. The Python scripts use the standard library. No package installation is performed automatically.

## 1. Recorded environment

The analyses used R 4.5.2 on macOS arm64 with the Accelerate BLAS and the RNG settings `Mersenne-Twister`, `Inversion`, `Rejection`.

| Package | Version | Additional pin |
|---|---|---|
| bkmr | 0.2.2.9000 | `jenfb/bkmr` commit `45413e338a316362d629f53bd2a917c4bf485c1e` |
| posterior | 1.7.0 | |
| tmvtnorm | 1.7 | |
| mvtnorm | 1.4-2 (R reports `1.4.2`) | |

The full recorded dependency environment is in `environment/`. A matching version is not a guarantee of matching compiled code, BLAS behavior, or floating-point trajectories. Use the recorded environment when exact matching is required. The runners intentionally stop if identity checks fail; investigate the differences instead of relaxing checks to claim reproduction. No claim is made that the complete study has been rerun on another operating system.

## 2. Quick checks of stored results

```sh
python3 Rev1/code/verify_results.py
```

This needs base R, but does not need bkmr, posterior, or any historical draws. Outputs are temporary and published CSVs are unchanged.

To keep an independent h aggregation table:

```sh
Rscript Rev1/code/reaggregate_h.R Rev1/_work/h_check
```

To keep the coauthor summaries:

```sh
python3 Rev1/coauthor_additions/reproduce_summaries.py --output-dir Rev1/_work/coauthor_check
```

## 3. Prepare an isolated workspace

```sh
python3 Rev1/code/manage.py prepare --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py check --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py audit --workdir Rev1/_work/rerun
```

`prepare` requires a new directory and refuses to replace an existing one. It copies the pinned analysis code and CSV reference inputs into the historical relative layout. `check` verifies key package versions, and regenerates all 12 extension datasets and 144 starting states against the exported exact specifications without MCMC. `audit` reconstructs the 430 primary classifications, overlap/acceptance tables, and failure-case input geometry without fitting.

The prepared workspace includes copies of the published initialization CSVs for immediate comparison-table/figure regeneration:

```sh
python3 Rev1/code/manage.py comparisons --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py figures --workdir Rev1/_work/rerun
```

The figure outputs are under `Rev1/_work/rerun/initialization_extension_20260914/figures/`. This is regeneration from stored summaries, not a new MCMC analysis. `init-summary` replaces the workspace summary copies only after validating all newly generated chain outcomes.

## 4. Pilot: one dataset

Task 201 is the n=40 anchor case. The following runs four 10,000-iteration baseline chains, then fits again with the same seeds and performs all-point h postprocessing:

```sh
python3 Rev1/code/manage.py beta --task 201 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py h --task 201 --workdir Rev1/_work/rerun
```

Both routes require the original six-parameter diagnostic values to match the distributed references. The h runner also requires the core draws to match the freshly regenerated baseline, verifies retained-iteration alignment, and checks the first conditional draw/moments. The full set of h diagnostics is written as CSV; only the nested five h histories plus core parameter histories are retained in each local checkpoint to reduce storage. These are local generated files and are not included in the distribution.

The public h route uses the final v2 runner for every dataset. Historically, six pilot datasets used the earlier runner; both historical sources are retained. The fitting model, SamplePred call, seeds, and diagnostic thresholds are unchanged.

## 5. Full beta and h analyses

These are explicitly long-running steps; they are not triggered by `prepare`, `check`, or `verify_results.py`.

```sh
python3 Rev1/code/manage.py beta --all --workers 8 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py beta-summary --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py h --all --workers 8 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py h-summary --workdir Rev1/_work/rerun
```

The public baseline stage covers tasks 51–250 (all four sample sizes) because those regenerated baselines are required for h and the extension. Section 4's beta performance/MCSE aggregation still uses only the original 149 completed n=80/160/320 cases. Task 174 is a numerical non-completion and is excluded from inferential summaries, not silently counted as diagnostic non-pass.

Expected paths within the workspace:

- `revision_20260911/beta_estimand_by_n_and_primary_status.csv`: published `tables/beta_performance.csv`.
- `revision_20260911/beta_mcse_summary.csv`: published `tables/beta_mcse.csv`.
- `revision_20260911/parameter_flag_overlap.csv`: published `tables/parameter_overlap.csv`.
- `revision_20260911/lambda_acceptance_diagnostics.csv`: published `tables/lambda_acceptance.csv`.
- `h_false_allpoints_20260912/`: point-, task-, and sample-size-level h tables.

Budget hours for full runs; runtime varies by hardware. Eight concurrent tasks can use substantial RAM. Reduce `--workers` as needed; the default is one. The batch runner stops dispatching after a failed validation and lets active jobs finish. It logs each task under `run_logs/`. A compatible completed checkpoint can be reused; changed code or environment invalidates the fingerprint. Use a new workspace after changing the execution code.

After the full h run, the baseline beta summaries by task and chain (199 datasets, 796 chains) can also be regenerated:

```sh
python3 Rev1/code/manage.py baseline-summary --workdir Rev1/_work/rerun
```

## 6. Initial-state and iteration-budget experiment

This can be run independently of the full h analysis. First regenerate only the twelve selected baselines, then prepare the extension inputs:

```sh
python3 Rev1/code/manage.py beta --init-cases --workers 8 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py init-inputs --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py init --all --workers 8 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py init-summary --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py figures --workdir Rev1/_work/rerun
```

After all twelve extension inputs are prepared, `initial-values` can regenerate the scalar/minimum/maximum initial-value table:

```sh
python3 Rev1/code/manage.py initial-values --workdir Rev1/_work/rerun
```

The 144 job IDs, case selections, sampler seeds and initialization seeds are in `initialization_extension/registry.csv`, `cases.csv`, and `plan.json`. Each chain is run toward 100,000 iterations; failed updates are excluded. Prefixes are evaluated only if all four chains completed them. The selected datasets and thresholds are not changed during reruns.

A single-chain example, after the task-201 baseline pilot:

```sh
python3 Rev1/code/manage.py init-inputs --task 201 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py init --job task_201_original_chain_1 --workdir Rev1/_work/rerun
```

The original-start prefix checks compare against baseline draws newly generated in this workspace, not against historical raw draws that were omitted from the distribution. Published baseline diagnostics remain a separate numerical check. `init-summary` requires all 144 chain outcomes and checks the twelve original-start 10k diagnostic sets against the historical summaries.

## 7. Numerical failure investigation and saved-draw example

After `audit` has generated task-174 inputs, reproduce the uninstrumented error before tracing it:

```sh
python3 Rev1/code/manage.py beta --task 174 --workdir Rev1/_work/rerun
python3 Rev1/code/manage.py failure-trace --workdir Rev1/_work/rerun
```

The trace is a separate refit that records the error-time state without changing sampler inputs or consuming extra RNG draws. Matrix checks are written to `revision_20260911/failure_case/matrix_checks.csv`. The stage requires the original error to have reproduced; it does not assume that every environment will reproduce the same numerical failure.

The historical task-90 ESS/MCSE example can be recalculated only after task 90's baseline and initialization inputs, and the task-174 trace, have been regenerated:

```sh
Rscript Rev1/coauthor_additions/check_saved_draws.R Rev1/_work/rerun
```

This command uses the supplied workspace path. It does not read the author's original local paths.

## 8. Preservation

Published outputs are in `Rev1/`; regenerated outputs are under the chosen workspace. Keep these separate. Raw RDS files created by a rerun are local artifacts and should remain excluded from Git. The study's distributed plan/provenance notes describe the historical analysis, including historical exact-draw checks; see `PROVENANCE.md` for what can be verified from this lightweight distribution.
