# Results directory

The completed pure BKMR outputs are grouped in:

- `pure_bkmr_20260616_startvalues/`

This subfolder contains the full 8-worker run completed on 2026-06-16 at
23:56:10 JST.

Key files:

- `pure_bkmr_20260616_startvalues/all_task_results.csv`
- `pure_bkmr_20260616_startvalues/task_table.csv`
- `pure_bkmr_20260616_startvalues/overall_runtime.csv`
- `pure_bkmr_20260616_startvalues/parameter_diagnostics_by_task/`

Diagnostic convention:

- `rhat` is computed with `posterior::rhat()`, which reports the maximum of
  rank-normalized split R-hat and rank-normalized folded-split R-hat.
- `ess_bulk` and `ess_tail` are computed with `posterior::ess_bulk()` and
  `posterior::ess_tail()`.
- The strict Vehtari-style diagnostic screen uses R-hat <= 1.01 and both
  bulk-ESS and tail-ESS >= 400.
- The lenient R-hat <= 1.05 screen is retained only as a sensitivity or
  screening reference, not as the recommended primary standard.
