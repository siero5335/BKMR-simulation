script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "_pure_bkmr_core.R"))
ensure_bkmr_available()

root <- project_root(script_file)
results_dir <- file.path(root, "results")
source_dir <- results_dir
out_dir <- file.path(results_dir, "parameter_diagnostics_by_task")
workers <- as.integer(Sys.getenv("BKMR_WORKERS", "8"))
if (is.na(workers) || workers < 1L) {
  workers <- 8L
}

old_parameter_file <- function(row) {
  sprintf(
    "%s__%s__n%d__varsel_%s__iter_%d.csv",
    row$scenario,
    row$dataset_seed,
    row$n,
    tolower(as.character(row$varsel)),
    row$iter
  )
}

task_parameter_file <- function(row) {
  paste(sprintf("task_%03d", row$task_index), old_parameter_file(row), sep = "__")
}

source_parameter_file <- function(row) {
  task_file <- task_parameter_file(row)
  task_path <- file.path(source_dir, row$bundle, "parameter_diagnostics", task_file)
  if (file.exists(task_path)) {
    return(task_file)
  }
  old_parameter_file(row)
}

new_parameter_file <- function(row, source_file) {
  if (grepl("^task_[0-9]{3}__", source_file)) {
    return(source_file)
  }
  paste(sprintf("task_%03d", row$task_index), source_file, sep = "__")
}

all_results <- read.csv(file.path(results_dir, "all_task_results.csv"), stringsAsFactors = FALSE)
ok <- all_results[all_results$status == "ok", , drop = FALSE]
ok$source_parameter_file <- vapply(seq_len(nrow(ok)), function(i) source_parameter_file(ok[i, ]), character(1))
ok$source_parameter_path <- file.path(source_dir, ok$bundle, "parameter_diagnostics", ok$source_parameter_file)
dup_key <- paste(ok$bundle, ok$source_parameter_file, sep = "/")
ok$duplicate_parameter_file <- ave(rep(TRUE, nrow(ok)), dup_key, FUN = length) > 1L

nondup <- ok[!ok$duplicate_parameter_file, , drop = FALSE]
for (i in seq_len(nrow(nondup))) {
  dest <- file.path(out_dir, nondup$bundle[i],
                    new_parameter_file(nondup[i, ], nondup$source_parameter_file[i]))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(nondup$source_parameter_path[i])) {
    stop("Missing source parameter diagnostics: ", nondup$source_parameter_path[i], call. = FALSE)
  }
  file.copy(nondup$source_parameter_path[i], dest, overwrite = TRUE)
}

dups <- ok[ok$duplicate_parameter_file, , drop = FALSE]
write_csv(
  dups[, intersect(c("bundle", "scenario", "dataset_seed", "n", "varsel", "iter", "warmup",
                     "row_id", "seed_id", "task_index", "chain_seeds",
                     "initial_value_seed", "initial_value_task_seed", "source_parameter_file"),
                   names(dups))],
  file.path(out_dir, "duplicate_parameter_tasks.csv")
)

run_duplicate <- function(row) {
  chain_seeds <- as.integer(strsplit(row$chain_seeds, ";", fixed = TRUE)[[1]])
  initial_value_seed <- if ("initial_value_seed" %in% names(row)) {
    as.integer(row$initial_value_seed)
  } else {
    fixed_initial_value_seed()
  }
  fit_row <- fit_task(
    scenario = row$scenario,
    dataset_seed = row$dataset_seed,
    n = row$n,
    hfun = row$hfun,
    beta_true = row$beta_true,
    iter = row$iter,
    warmup = row$warmup,
    varsel = row$varsel,
    chain_seeds = chain_seeds,
    parameter_dir = file.path(out_dir, row$bundle),
    parameter_id = sprintf("task_%03d", row$task_index),
    initial_value_seed = initial_value_seed
  )
  fit_row$bundle <- row$bundle
  fit_row$row_id <- row$row_id
  fit_row$seed_id <- row$seed_id
  fit_row$task_index <- row$task_index
  fit_row
}

duplicate_rows <- split(dups, seq_len(nrow(dups)))
rerun_rows <- parallel::mclapply(duplicate_rows, run_duplicate, mc.cores = workers)
rerun_results <- do.call(rbind, rerun_rows)
write_csv(rerun_results, file.path(out_dir, "rerun_duplicate_task_results.csv"))

by_task_files <- list.files(out_dir, pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
by_task_files <- by_task_files[!basename(by_task_files) %in% c(
  "duplicate_parameter_tasks.csv",
  "rerun_duplicate_task_results.csv"
)]
write_csv(
  data.frame(
    all_ok_tasks = nrow(ok),
    copied_nonduplicate_tasks = nrow(nondup),
    rerun_duplicate_tasks = nrow(dups),
    by_task_parameter_files = length(by_task_files),
    workers = workers,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "by_task_parameter_diagnostics_summary.csv")
)

message("Created ", length(by_task_files), " by-task parameter diagnostics files.")
