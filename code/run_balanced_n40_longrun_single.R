source(file.path(dirname(script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "_pure_bkmr_core.R"))
ensure_bkmr_available()
root <- project_root(script_path)
out_dir <- file.path(root, "results", "balanced_n40_longrun_single")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chain_seeds <- fixed_chain_seeds()
initial_value_seed <- fixed_initial_value_seed()
rows <- lapply(c(50000L, 100000L), function(iter) {
  fit_task("balanced_bkmr_n40_single_longrun", 303L, 40L, 2L, 0.5,
           iter, iter / 2L, FALSE, chain_seeds,
           parameter_dir = file.path(out_dir, "parameter_diagnostics"),
           initial_value_seed = initial_value_seed)
})
details <- do.call(rbind, rows)
write_csv(details, file.path(out_dir, "summary.csv"))
write_csv(data.frame(dataset_seed = 303L, n = 40L,
                     chain_seeds = paste(chain_seeds, collapse = ";"),
                     initial_value_seed = initial_value_seed,
                     starting_values = bkmr_starting_values_description()),
          file.path(out_dir, "meta.csv"))
