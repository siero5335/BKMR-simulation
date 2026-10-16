# Packaging validation — 16 October 2026

Checks performed on the author's macOS arm64 / R 4.5.2 environment:

- All pre-existing result CSVs are byte-for-byte unchanged from the final revision folder.
- Every historical source file matches its SHA-256 entry in `code/SOURCE_INDEX.csv`.
- All distributed R source/input text files parse and all Python scripts pass syntax inspection.
- A fresh isolated workspace was prepared without any historical draws in its inputs.
- All twelve simulated datasets and all 144 initial states regenerated exactly against the exported historical input specifications.
- The original 431/430/1 task counts, 430 primary classifications, overlap counts, and failure-case input geometry checks passed.
- The four published initialization comparison tables and the coauthor summary JSON/tables regenerated from stored summaries.
- All 29,680 h point records were independently reaggregated; the primary mapping and nested-five/all-point comparison agreed.
- Task 201 (n=40) was refitted with four 10,000-iteration chains. Its original six-parameter diagnostic check passed. This took about 25 seconds on the test machine.
- Task 201 was then refitted through the public all-point h route. All 40 h diagnostic rows matched the historical point-level records within 1e-10; core draws matched the fresh baseline exactly. This took about 26 seconds.
- One extension job, `task_201_original_chain_1`, completed 100,000 iterations through the public runner. Its baseline prefix check passed. A packaging-only comparison against the author's locally retained historical trajectory also found an exact match across all 100,000 iterations. That historical trajectory is not included in the package.
- The beta summary adapter was checked using the existing historical 149 completed-task checkpoints linked into the test workspace. The beta performance and MCSE tables matched. This was a reaggregation test, not 149 fresh MCMC refits; none of those checkpoints or links is included in the package.
- All four files in `tables/` matched their regenerated counterparts. Initialization figures regenerated from published summary CSVs.

The full 200-task h rerun and all 144 extension jobs were not rerun during packaging. The public runtime's full h/extension validation stages therefore remain untested as complete new batches. This package does not claim independent reimplementation of posterior's R-hat/ESS algorithms or verification on a different operating system. The historical analysis/validation and the new public rerun checks are distinguished in `PROVENANCE.md`.

No raw draw files, models, RDS checkpoints, test workspaces, logs, or symlinks are distributed. No git commit or push was performed during preparation of the downloadable package.
