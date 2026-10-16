# Provenance and public rerun adaptations

The published numeric CSVs were copied from the author's final revision materials without changing their values. The extra all-observed-point h CSV and associated core diagnostic CSV were taken from the completed September all-point analysis and its verification addendum. The twelve initialization specifications were exported from the historical inputs using an explicit allowlist: task metadata, simulated data, starting-state sets, baseline diagnostic summaries, sampler seeds, initialization seed, and RNG kind. No reference draws or chain histories were exported.

The scripts in `code/original/` are byte-for-byte copies of the historical analytical scripts. `code/SOURCE_INDEX.csv` gives their source-relative locations and hashes. Other historical orchestration scripts, draft-document builders unrelated to numerical aggregation, and logs are not required by the public entry point. `make_documents.py` is retained in the historical source because its first section produced the published initialization comparison tables; the public runtime executes only that numerical section, not its manuscript/archive-writing sections.

## Public adapters

| Component | Adaptation | Preserved checks |
|---|---|---|
| Workspace preparation | Rebuild historical relative paths in a new directory from distributed code and CSVs; split consolidated diagnostics into per-task references | Refuse existing target directories; published inputs remain separate |
| Original-result audit | Omit the scan of the author's unrelated `old/` RDS files and the resulting local inventory narrative | Original 431/430/1 task counts, 430 primary labels, diagnostic and geometry calculations |
| Baseline refit | Also allow the n=40 balanced-ablation subset; run all four sample sizes as needed by h/extension | Original data/start generators, seeds, model, package pin, diagnostics, fingerprint, CSV-level identity gate |
| All-point h | Use a freshly generated task-specific baseline for the package fingerprint and exact core-draw comparison; use final v2 runner for all tasks | Original fits, postprocessing seeds, SamplePred, independent moments, first-draw replay, iteration alignment, original diagnostic comparison |
| h validation | Compare against the new baseline in the same workspace; omit the old local source-manifest command | All chain/core/h checks and aggregate expected counts |
| Extension input preparation | Generate baseline reference draws locally; supply exact data/starts and historical baseline diagnostics from text specifications | Exact input/start checks; baseline identity gate |
| Extension fitting | Historical `run_chain.R` copied unchanged to runtime | Job/input/plan/package fingerprints, original-start prefix matching to regenerated baseline, complete-iteration recovery |
| Extension summary | Omit the historical source-manifest command; label prefix checks as regenerated-baseline checks | 144 job outcomes, fingerprints, starts, prefix checks, twelve historical baseline diagnostic comparisons |
| Comparison tables | Extract only the numerical section from historical `make_documents.py` | Same aggregation formulas and input CSVs |
| Public orchestration | Configurable workers/Rscript; no draft-document generation or implicit launch of all MCMC | Explicit stages, task logs, stop dispatch on validation failure |
| Saved-draw example | Accept a user-specified rerun workspace | Same numerical formulas; raw inputs must have been regenerated |

The historical source and runtime adapters are distinguishable in the directory names. The public adapter files are not claimed to have been frozen before the historical experiments. Historical selection/freeze information is in `initialization_extension/plan.json` and `PROVENANCE_JA.md`.

## Exact trajectory checks

The historical full analysis used saved previous draws to check exact reproduction. Those draws remain outside this public package. Public reruns first check regenerated baseline diagnostics against the published original references, then compare the h/extension trajectories with those freshly generated baselines. These are two different checks. A public rerun cannot directly prove identity to omitted historical trajectories merely by passing a comparison between two new runs.

Floating-point or dependency differences can cause an identity check to fail. Failures are surfaced; results are not silently relabelled as reproductions. The full historical h analysis included six earlier pilot executions; the new uniform v2 route has been tested on the n=40 anchor, not rerun for all 200 datasets during packaging.
