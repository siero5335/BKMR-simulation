# BKMR simulation and MCMC diagnostics

Code and numerical results for the original BKMR simulation study and the additional analyses undertaken during the first manuscript revision.

| Material | Location |
|---|---|
| Original simulation code | [`code/`](code/) |
| Original simulation results | [`results/`](results/) |
| Revision analyses, lightweight results, and reproduction instructions | [`Rev1/`](Rev1/) |

The original design contains 431 tasks: 430 completed fits and one numerical non-completion (`balanced_bkmr_n320_varsel_false`, dataset seed 1825111885; task 174). Task-specific data, initialization, and sampler seeds are recorded in the task table; do not assume that a single set of seeds applies to every study bundle.

The primary diagnostic criterion requires all six monitored parameters (beta, lambda, and four r parameters) to meet R-hat <= 1.01, bulk-ESS >= 400, and tail-ESS >= 400. The revision analyses retain this criterion and the original task classifications.

## First revision

`Rev1/` contains additional beta estimand/MCSE summaries, all-observed-point h diagnostics, initialization and iteration-budget comparisons, numerical failure investigation, and supporting numerical checks. The README maps these materials to the Supplementary Information.

For a quick numerical check without MCMC, from the repository root:

```sh
python3 Rev1/code/verify_results.py
```

This requires Python 3.9+ and base R. Full refits require the recorded R packages and are explicitly invoked through the isolated-workspace runner described in [`Rev1/REPRODUCING.md`](Rev1/REPRODUCING.md).

No raw MCMC trajectories, fitted model objects, or manuscripts are included in the revision directory. The small simulated-data and initialization specifications are included; reruns generate their own trajectories locally. Published summary files remain separate from regenerated outputs.
