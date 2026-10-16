Numerical additions from coauthor review (13 October 2026), supporting Supplemental Table 5.2, descriptive endpoint-MCSE counts in SI 4.2, event-prevalence comparisons, and task-specific sampler-seed clarification.

From the repository root:
python3 Rev1/coauthor_additions/reproduce_summaries.py --output-dir Rev1/_work/coauthor_check

Only Python's standard library and distributed CSVs are needed. This regenerates numeric_checks.json, reference_comparisons.csv (54 evaluations), and chain_ranges.csv (144 evaluations) in a separate output directory. It does not overwrite the published reference tables or fit MCMC. Chain means are independently joined to fit summaries. Reference differences are absolute differences divided by the evaluated fit's posterior SD. Evaluations are selected and dependent.

beta_estimand_by_task.csv, task_table.csv, and all_task_results.csv are unchanged existing study outputs. To verify all reproduced values against the distributed references, use python3 Rev1/code/verify_results.py.

check_saved_draws.R and saved_draw_checks.txt document a separate task-90/task-174 check. The text is the historical check output. To rerun the script, first regenerate the required task-90 baseline/extension inputs and task-174 failure trace as described in ../REPRODUCING.md, then supply the new workspace path:
Rscript Rev1/coauthor_additions/check_saved_draws.R Rev1/_work/rerun

No historical RDS or raw draws are distributed. This R check uses posterior and already regenerated saved draws; it does not itself rerun MCMC.
