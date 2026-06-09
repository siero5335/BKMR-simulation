# Public bkmr Verification Memo

This directory contains an independent verification subset for the currently
installed public GitHub `bkmr` package.

## Environment identity

- Installed package path:
  - `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bkmr`
- Installed version:
  - `0.2.2.9000`
- Recorded `RemoteSha`:
  - `45413e338a316362d629f53bd2a917c4bf485c1e`

See `verification_env.txt` for the full environment record and the installed
`NEWS.md` excerpt.

## Verification steps included

1. Main benchmark, exact upstream only
   - `n = 40`, `M = 4`, `iter = 10000`, `warmup = 5000`, `chains = 4`,
     `varsel = TRUE`
2. `n = 80` crash check
   - `iter = 10000`, `varsel = FALSE`
3. Seed robustness
   - three additional independent chain-seed sets for Step 1
4. Long-chain sanity check
   - `n = 40`, `iter = 100000`, `warmup = 50000`, `varsel = FALSE`
5. Extended seed robustness
   - 10 independent seed sets for the Step 1 configuration
6. Sample-size sweep
   - exact-only, `varsel = FALSE`, `n = 40, 80, 160, 320`
7. Exact long-run check
   - exact-only, `varsel = FALSE`, `50k/25k` and `100k/50k`
8. Exact long-run seed robustness
   - exact-only, `varsel = FALSE`, `50k/25k` and `100k/50k`, 3 seed sets

## Main numerical result

The public installed package is qualitatively consistent with the paper's
reported concerns.

- Step 1 main exact run:
  - `max_rhat = 2.747`
  - `min_ess = 40.56`
  - hard convergence failed
- Step 3 robustness runs:
  - all three additional seed sets also failed hard convergence
  - `max_rhat` ranged from about `1.59` to `3.05`
  - `min_ess` ranged from about `30.2` to `59.0`
- Step 5 extended seed check:
  - all 10 seed sets completed without crash
  - all 10 seed sets still failed hard convergence
  - this rules out a single unlucky seed explanation
- Step 6 sample-size sweep:
  - `n = 40` gave `max_rhat = 2.662`, `min_ess = 54.47`
  - `n = 80, 160, 320` each returned
    `missing value where TRUE/FALSE needed`
- Step 4 long-chain exact run:
  - `max_rhat = 1.283`
  - `min_ess = 89.20`
  - hard convergence still failed
- Step 7 exact long-run extension:
  - `50k/25k`: `max_rhat = 1.513`, `min_ess = 52.85`
  - `100k/50k`: `max_rhat = 1.283`, `min_ess = 89.20`
  - both still failed hard convergence
- Step 8 exact long-run seed robustness:
  - `50k/25k`, 3 seed sets:
    - `crash_rate = 0`
    - `hard_failure_rate = 1`
    - `median_max_rhat = 1.611`
    - `worst_max_rhat = 1.734`
    - `median_min_ess = 57.77`
  - `100k/50k`, 3 seed sets:
    - `crash_rate = 0`
    - `hard_failure_rate = 1`
    - `median_max_rhat = 1.160`
    - `worst_max_rhat = 1.180`
    - `median_min_ess = 77.36`
- Step 2 crash check:
  - `n = 80`, `varsel = FALSE` returned
    `missing value where TRUE/FALSE needed`

## Files

- `verification_env.txt`
  - package provenance, installed `NEWS.md` excerpt, session info
- `main_verification_summary.csv`
  - compact table for Steps 1, 3, and 4
- `crash_verification.txt`
  - exact result for Step 2
- `seed_table.csv`
  - chain seeds used in Steps 1, 3, and 4
- `seed10_summary.csv`
  - compact 10-seed exact-only robustness summary
- `seed10_run_details.csv`
  - per-seed diagnostics for the 10-seed check
- `seed10_seed_table.csv`
  - chain seeds used in the 10-seed check
- `n_sweep_summary.csv`
  - public exact-only `n = 40, 80, 160, 320` sweep
- `n_sweep_meta.csv`
  - settings used for the public exact-only `n` sweep
- `longrun_50k_100k_summary.csv`
  - public exact-only `50k/25k` and `100k/50k` comparison
- `longrun_50k_100k_meta.csv`
  - settings used for the public exact-only long-run comparison
- `longrun_seed3_summary.csv`
  - 3-seed robustness summary for the public exact-only `50k/25k` and
    `100k/50k` runs
- `longrun_seed3_run_details.csv`
  - per-seed diagnostics for the public exact-only longrun robustness check
- `longrun_seed3_seed_table.csv`
  - chain seeds used in the public exact-only longrun robustness check
- `longrun_seed3_meta.csv`
  - settings used for the public exact-only longrun robustness check
- `meta.csv`
  - global verification settings
- `verify_on_public_bkmr.R`
  - standalone reproduction script
- `verify_public_bkmr_seed10.R`
  - standalone script for the 10-seed robustness extension
- `verify_public_bkmr_n_sweep.R`
  - standalone script for the public exact-only `n` sweep
- `verify_public_bkmr_longrun_50k_100k.R`
  - standalone script for the public exact-only `50k/25k` and `100k/50k` check
- `verify_public_bkmr_longrun_seed3.R`
  - standalone script for the public exact-only longrun 3-seed robustness check
