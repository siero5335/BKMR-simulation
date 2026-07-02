source(file.path(dirname(script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "_pure_bkmr_core.R"))
ensure_bkmr_available()
root <- project_root(script_path)
out_dir <- file.path(root, "results", "balanced_longrun_seed3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chain_seeds <- fixed_chain_seeds()
initial_value_seed <- fixed_initial_value_seed()
seeds <- dataset_seeds(20260522L, 3L)
grid <- expand.grid(
  n = c(40L, 80L),
  iter = c(50000L, 100000L),
  KEEP.OUT.ATTRS = FALSE
)
rows <- list()
for (i in seq_len(nrow(grid))) {
  for (seed in seeds) {
    rows[[length(rows) + 1L]] <- fit_task(
      scenario = sprintf("balanced_bkmr_n%d_varsel_false_longrun", grid$n[i]),
      dataset_seed = seed,
      n = grid$n[i],
      hfun = 2L,
      beta_true = 0.5,
      iter = grid$iter[i],
      warmup = grid$iter[i] / 2L,
      varsel = FALSE,
      chain_seeds = chain_seeds,
      parameter_dir = file.path(out_dir, "parameter_diagnostics"),
      initial_value_seed = initial_value_seed
    )
  }
}
details <- do.call(rbind, rows)
summary <- summarize_sweep(details, c("scenario", "n", "varsel", "hfun", "beta_true", "iter", "warmup"))
write_csv(details, file.path(out_dir, "run_details.csv"))
write_csv(summary, file.path(out_dir, "summary.csv"))
write_csv(data.frame(dataset_seed_source = 20260522L, K = 3L,
                     chain_seeds = paste(chain_seeds, collapse = ";"),
                     initial_value_seed = initial_value_seed,
                     starting_values = bkmr_starting_values_description()),
          file.path(out_dir, "meta.csv"))
