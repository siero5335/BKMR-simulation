script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "_pure_bkmr_core.R"))
ensure_bkmr_available()

root <- project_root(script_file)
results_dir <- file.path(root, "results")
run_started_at <- Sys.time()
workers <- as.integer(Sys.getenv("BKMR_WORKERS", "8"))
if (is.na(workers) || workers < 1L) {
  workers <- 8L
}

task_row <- function(bundle, scenario, dataset_seed, n, hfun, beta_true, iter, warmup,
                     varsel, chain_seeds, row_id = NA_character_, seed_id = NA_integer_) {
  data.frame(
    bundle = bundle,
    scenario = scenario,
    dataset_seed = dataset_seed,
    n = n,
    hfun = hfun,
    beta_true = beta_true,
    iter = iter,
    warmup = warmup,
    varsel = varsel,
    chain_seeds = paste(chain_seeds, collapse = ";"),
    row_id = row_id,
    seed_id = seed_id,
    stringsAsFactors = FALSE
  )
}

fixed <- fixed_chain_seeds()
tasks <- list()

balanced_seeds <- dataset_seeds(20260522L, 50L)
balanced_grid <- data.frame(
  scenario = c("balanced_bkmr_n40_varsel_true", "balanced_bkmr_n80_varsel_false",
               "balanced_bkmr_n160_varsel_false", "balanced_bkmr_n320_varsel_false"),
  n = c(40L, 80L, 160L, 320L),
  varsel = c(TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
for (g in seq_len(nrow(balanced_grid))) {
  for (seed in balanced_seeds) {
    tasks[[length(tasks) + 1L]] <- task_row(
      "balanced_main", balanced_grid$scenario[g], seed, balanced_grid$n[g],
      2L, 0.5, 10000L, 5000L, balanced_grid$varsel[g], fixed
    )
  }
}

for (seed in balanced_seeds) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "balanced_ablation", "balanced_bkmr_n40_varsel_false", seed, 40L,
    2L, 0.5, 10000L, 5000L, FALSE, fixed
  )
}

for (seed in dataset_seeds(20260522L, 3L)) {
  for (n in c(40L, 80L)) {
    for (iter in c(50000L, 100000L)) {
      tasks[[length(tasks) + 1L]] <- task_row(
        "balanced_longrun_seed3",
        sprintf("balanced_bkmr_n%d_varsel_false_longrun", n),
        seed, n, 2L, 0.5, iter, iter / 2L, FALSE, fixed
      )
    }
  }
}

for (iter in c(50000L, 100000L)) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "balanced_n40_longrun_single", "balanced_bkmr_n40_single_longrun",
    303L, 40L, 2L, 0.5, iter, iter / 2L, FALSE, fixed
  )
}

for (seed in dataset_seeds(20260522L, 20L)) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "default_support", "default_bkmr_n40_varsel_true", seed, 40L,
    3L, 2, 10000L, 5000L, TRUE, fixed
  )
  tasks[[length(tasks) + 1L]] <- task_row(
    "default_support", "default_bkmr_n80_varsel_false", seed, 80L,
    3L, 2, 10000L, 5000L, FALSE, fixed
  )
}

dat40_for_seed <- function(seed) {
  make_bkmr_dataset(seed, 40L, 3L, 2)
}
set.seed(303L)
invisible(bkmr::SimData(n = 40, M = 4, family = "binomial", hfun = 3L, beta.true = 2, ind = 1:2))
public_main_chain_seeds <- sample.int(.Machine$integer.max, 4L)
set.seed(20260420L)
public_seed_sets_3 <- replicate(3L, sample.int(.Machine$integer.max, 4L))
set.seed(20260420L)
public_seed_sets_10 <- replicate(10L, sample.int(.Machine$integer.max, 4L))

tasks[[length(tasks) + 1L]] <- task_row(
  "public_main", "step1_main_exact_varsel_true", 303L, 40L,
  3L, 2, 10000L, 5000L, TRUE, public_main_chain_seeds, "step1_main"
)
for (i in seq_len(ncol(public_seed_sets_3))) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_main", sprintf("step3_seed_%d_exact_varsel_true", i), 303L, 40L,
    3L, 2, 10000L, 5000L, TRUE, public_seed_sets_3[, i],
    sprintf("step3_seed_%d", i), i
  )
}
tasks[[length(tasks) + 1L]] <- task_row(
  "public_main", "step4_long_exact_varsel_false", 303L, 40L,
  3L, 2, 100000L, 50000L, FALSE, public_main_chain_seeds, "step4_long"
)

for (i in seq_len(ncol(public_seed_sets_10))) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_seed10", "public_bkmr_seed10_varsel_true", 303L, 40L,
    3L, 2, 10000L, 5000L, TRUE, public_seed_sets_10[, i],
    sprintf("seed10_%02d", i), i
  )
}

for (seed in dataset_seeds(20260514L, 50L)) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_dataset_sweep", "public_bkmr_dataset_sweep_n40_varsel_true",
    seed, 40L, 3L, 2, 10000L, 5000L, TRUE, fixed
  )
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_dataset_sweep", "public_bkmr_dataset_sweep_n80_varsel_false",
    seed, 80L, 3L, 2, 10000L, 5000L, FALSE, fixed
  )
}

for (n in c(40L, 80L, 160L, 320L)) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_n_sweep", "public_bkmr_n_sweep_varsel_false",
    303L, n, 3L, 2, 10000L, 5000L, FALSE, fixed
  )
}

for (iter in c(50000L, 100000L)) {
  tasks[[length(tasks) + 1L]] <- task_row(
    "public_longrun_50k_100k", "public_bkmr_longrun_varsel_false",
    303L, 40L, 3L, 2, iter, iter / 2L, FALSE, fixed
  )
}

for (i in seq_len(ncol(public_seed_sets_3))) {
  for (iter in c(50000L, 100000L)) {
    tasks[[length(tasks) + 1L]] <- task_row(
      "public_longrun_seed3", "public_bkmr_longrun_seed3_varsel_false",
      303L, 40L, 3L, 2, iter, iter / 2L, FALSE, public_seed_sets_3[, i],
      sprintf("seed_%d_iter_%d", i, iter), i
    )
  }
}

task_table <- do.call(rbind, tasks)
task_table$task_index <- seq_len(nrow(task_table))
write_csv(task_table, file.path(results_dir, "task_table.csv"))

message("Running ", nrow(task_table), " tasks with ", workers, " workers.")

run_task <- function(task) {
  chain_seeds <- as.integer(strsplit(task$chain_seeds, ";", fixed = TRUE)[[1]])
  parameter_dir <- file.path(results_dir, task$bundle, "parameter_diagnostics")
  row <- fit_task(
    scenario = task$scenario,
    dataset_seed = task$dataset_seed,
    n = task$n,
    hfun = task$hfun,
    beta_true = task$beta_true,
    iter = task$iter,
    warmup = task$warmup,
    varsel = task$varsel,
    chain_seeds = chain_seeds,
    parameter_dir = parameter_dir,
    parameter_id = sprintf("task_%03d", task$task_index)
  )
  row$bundle <- task$bundle
  row$row_id <- task$row_id
  row$seed_id <- task$seed_id
  row$task_index <- task$task_index
  row
}

task_list <- split(task_table, seq_len(nrow(task_table)))
rows <- parallel::mclapply(task_list, run_task, mc.cores = workers)
all_results <- do.call(rbind, rows)
write_csv(all_results, file.path(results_dir, "all_task_results.csv"))

errors <- all_results[all_results$status != "ok", c(
  "bundle", "scenario", "dataset_seed", "n", "varsel", "iter", "warmup",
  "label_sum", "error_message"
)]
write_csv(errors, file.path(results_dir, "pure_bkmr_error_tasks.csv"))

write_group <- function(bundle, file, data) {
  write_csv(data, file.path(results_dir, bundle, file))
}

balanced_main <- all_results[all_results$bundle == "balanced_main", , drop = FALSE]
balanced_main$dataset_sweep_seed <- 20260522L
write_group("balanced_main", "run_details.csv", balanced_main)
write_group("balanced_main", "summary.csv",
            summarize_sweep(balanced_main, c("scenario", "n", "varsel", "dataset_sweep_seed")))
write_group("balanced_main", "meta.csv",
            data.frame(dataset_sweep_seed = 20260522L, K = 50L,
                       iter = 10000L, warmup = 5000L,
                       chain_seeds = paste(fixed, collapse = ";")))

balanced_ablation <- all_results[all_results$bundle == "balanced_ablation", , drop = FALSE]
balanced_ablation$dataset_sweep_seed <- 20260522L
write_group("balanced_ablation", "run_details.csv", balanced_ablation)
write_group("balanced_ablation", "summary.csv",
            summarize_sweep(balanced_ablation, c("scenario", "n", "varsel", "dataset_sweep_seed")))
write_group("balanced_ablation", "meta.csv",
            data.frame(dataset_sweep_seed = 20260522L, K = 50L,
                       iter = 10000L, warmup = 5000L,
                       chain_seeds = paste(fixed, collapse = ";")))

balanced_longrun <- all_results[all_results$bundle == "balanced_longrun_seed3", , drop = FALSE]
write_group("balanced_longrun_seed3", "run_details.csv", balanced_longrun)
write_group("balanced_longrun_seed3", "summary.csv",
            summarize_sweep(balanced_longrun, c("scenario", "n", "varsel", "hfun", "beta_true", "iter", "warmup")))
write_group("balanced_longrun_seed3", "meta.csv",
            data.frame(dataset_seed_source = 20260522L, K = 3L,
                       chain_seeds = paste(fixed, collapse = ";")))

balanced_single <- all_results[all_results$bundle == "balanced_n40_longrun_single", , drop = FALSE]
write_group("balanced_n40_longrun_single", "summary.csv", balanced_single)
write_group("balanced_n40_longrun_single", "meta.csv",
            data.frame(dataset_seed = 303L, n = 40L,
                       chain_seeds = paste(fixed, collapse = ";")))

default_support <- all_results[all_results$bundle == "default_support", , drop = FALSE]
default_support$dataset_sweep_seed <- 20260522L
write_group("default_support", "run_details.csv", default_support)
write_group("default_support", "summary.csv",
            summarize_sweep(default_support, c("scenario", "n", "varsel", "dataset_sweep_seed")))
write_group("default_support", "meta.csv",
            data.frame(dataset_sweep_seed = 20260522L, K = 20L,
                       iter = 10000L, warmup = 5000L,
                       chain_seeds = paste(fixed, collapse = ";")))

public_main <- all_results[all_results$bundle == "public_main", , drop = FALSE]
public_main$dataset_seed <- 303L
public_main$public_bkmr_version <- as.character(utils::packageVersion("bkmr"))
public_main$public_bkmr_sha <- utils::packageDescription("bkmr")$RemoteSha
write_group("public_verification", "main_verification_summary.csv", public_main)

public_seed10 <- all_results[all_results$bundle == "public_seed10", , drop = FALSE]
ok <- public_seed10$status == "ok"
completed <- public_seed10[ok, , drop = FALSE]
write_group("public_verification", "seed10_run_details.csv", public_seed10)
write_group("public_verification", "seed10_summary.csv", data.frame(
  scenario = "public_bkmr_seed10_varsel_true",
  model = "exact_public",
  dataset_seed = 303L,
  seed_sweep_seed = 20260420L,
  n_seed_sets = 10L,
  n_completed = sum(ok),
  crash_rate = mean(!ok),
  hard_failure_rate = if (sum(ok) == 0L) NA_real_ else mean(!completed$converged_hard),
  tuning_issue_rate = if (sum(ok) == 0L) NA_real_ else mean(completed$tuning_issue),
  median_max_rhat = if (sum(ok) == 0L) NA_real_ else stats::median(completed$max_rhat, na.rm = TRUE),
  worst_max_rhat = if (sum(ok) == 0L) NA_real_ else max(completed$max_rhat, na.rm = TRUE),
  median_min_ess = if (sum(ok) == 0L) NA_real_ else stats::median(completed$min_ess, na.rm = TRUE),
  worst_min_ess = if (sum(ok) == 0L) NA_real_ else min(completed$min_ess, na.rm = TRUE),
  median_min_bulk_ess = if (sum(ok) == 0L) NA_real_ else stats::median(completed$min_bulk_ess, na.rm = TRUE),
  worst_min_bulk_ess = if (sum(ok) == 0L) NA_real_ else min(completed$min_bulk_ess, na.rm = TRUE),
  median_min_tail_ess = if (sum(ok) == 0L) NA_real_ else stats::median(completed$min_tail_ess, na.rm = TRUE),
  worst_min_tail_ess = if (sum(ok) == 0L) NA_real_ else min(completed$min_tail_ess, na.rm = TRUE),
  stringsAsFactors = FALSE
))

public_dataset <- all_results[all_results$bundle == "public_dataset_sweep", , drop = FALSE]
public_dataset$dataset_sweep_seed <- 20260514L
write_group("public_verification", "dataset_sweep_run_details.csv", public_dataset)
write_group("public_verification", "dataset_sweep_summary.csv",
            summarize_sweep(public_dataset, c("scenario", "dataset_sweep_seed")))

public_n <- all_results[all_results$bundle == "public_n_sweep", , drop = FALSE]
write_group("public_verification", "n_sweep_summary.csv", public_n)

public_long <- all_results[all_results$bundle == "public_longrun_50k_100k", , drop = FALSE]
write_group("public_verification", "longrun_50k_100k_summary.csv", public_long)

public_long3 <- all_results[all_results$bundle == "public_longrun_seed3", , drop = FALSE]
write_group("public_verification", "longrun_seed3_run_details.csv", public_long3)
write_group("public_verification", "longrun_seed3_summary.csv",
            summarize_sweep(public_long3, c("scenario", "iter", "warmup")))

env_path <- file.path(results_dir, "public_verification", "verification_env.txt")
sink(env_path)
cat("Public bkmr verification environment\n")
cat("Package path:", find.package("bkmr"), "\n")
cat("Version:", as.character(utils::packageVersion("bkmr")), "\n")
cat("RemoteSha:", utils::packageDescription("bkmr")$RemoteSha, "\n\n")
cat("posterior package path:", find.package("posterior"), "\n")
cat("posterior version:", as.character(utils::packageVersion("posterior")), "\n\n")
print(sessionInfo())
sink()

crash_check <- tryCatch({
  dat80 <- make_bkmr_dataset(303L, 80L, 3L, 2)
  set.seed(public_main_chain_seeds[1])
  invisible(utils::capture.output({
    bkmr::kmbayes(y = dat80$y, Z = dat80$Z, X = dat80$X, iter = 10000L,
                  family = "binomial", varsel = FALSE, verbose = FALSE)
  }))
  list(status = "ran_to_completion", message = NA_character_)
}, error = function(e) {
  list(status = "error", message = conditionMessage(e))
})
writeLines(c(
  "n = 80 numerical crash check",
  paste("dataset_seed:", 303L),
  paste("chain_seed:", public_main_chain_seeds[1]),
  paste("iter:", 10000L),
  paste("warmup:", 5000L),
  paste("status:", crash_check$status),
  paste("message:", crash_check$message)
), file.path(results_dir, "public_verification", "crash_verification.txt"))

run_finished_at <- Sys.time()
write_csv(data.frame(
  started_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
  finished_at = format(run_finished_at, "%Y-%m-%d %H:%M:%S %Z"),
  elapsed_sec = as.numeric(difftime(run_finished_at, run_started_at, units = "secs")),
  workers = workers,
  n_tasks = nrow(task_table),
  n_ok = sum(all_results$status == "ok"),
  n_error = sum(all_results$status != "ok"),
  diagnostic_package = "posterior",
  diagnostic_package_version = as.character(utils::packageVersion("posterior")),
  rhat_definition = "maximum of rank-normalized split R-hat and rank-normalized folded-split R-hat",
  rhat_threshold = 1.01,
  ess_threshold = 400,
  starting_values = "explicit binomial starts: h.hat=1, beta=0, sigsq.eps=1, r=1, lambda=10, delta=1, ystar=+/-0.5",
  stringsAsFactors = FALSE
), file.path(results_dir, "overall_runtime.csv"))

message("Completed tasks: ", sum(all_results$status == "ok"), " ok / ",
        sum(all_results$status != "ok"), " error.")
