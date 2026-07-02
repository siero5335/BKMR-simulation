source(file.path(dirname(script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "_pure_bkmr_core.R"))
ensure_bkmr_available()
root <- project_root(script_path)
out_dir <- file.path(root, "results", "balanced_main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

K <- 50L
iter <- 10000L
warmup <- 5000L
sweep_seed <- 20260522L
chain_seeds <- fixed_chain_seeds()
initial_value_seed <- fixed_initial_value_seed()
seeds <- dataset_seeds(sweep_seed, K)
grid <- data.frame(
  scenario = c("balanced_bkmr_n40_varsel_true", "balanced_bkmr_n80_varsel_false",
               "balanced_bkmr_n160_varsel_false", "balanced_bkmr_n320_varsel_false"),
  n = c(40L, 80L, 160L, 320L),
  varsel = c(TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

rows <- list()
for (g in seq_len(nrow(grid))) {
  for (seed in seeds) {
    rows[[length(rows) + 1L]] <- fit_task(
      scenario = grid$scenario[g],
      dataset_seed = seed,
      n = grid$n[g],
      hfun = 2L,
      beta_true = 0.5,
      iter = iter,
      warmup = warmup,
      varsel = grid$varsel[g],
      chain_seeds = chain_seeds,
      parameter_dir = file.path(out_dir, "parameter_diagnostics"),
      initial_value_seed = initial_value_seed
    )
  }
}
details <- do.call(rbind, rows)
details$dataset_sweep_seed <- sweep_seed
summary <- summarize_sweep(details, c("scenario", "n", "varsel", "dataset_sweep_seed"))
write_csv(details, file.path(out_dir, "run_details.csv"))
write_csv(summary, file.path(out_dir, "summary.csv"))
write_csv(data.frame(dataset_sweep_seed = sweep_seed, K = K, iter = iter,
                     warmup = warmup, chain_seeds = paste(chain_seeds, collapse = ";"),
                     initial_value_seed = initial_value_seed,
                     starting_values = bkmr_starting_values_description()),
          file.path(out_dir, "meta.csv"))
