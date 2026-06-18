source(file.path(dirname(script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "_pure_bkmr_core.R"))
ensure_bkmr_available()
root <- project_root(script_path)
out_dir <- file.path(root, "results", "public_verification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_existing_data <- function(dat, scenario, row_id, iter, warmup, varsel, chain_seeds) {
  fit_obj <- NULL
  elapsed <- system.time({
    fit_obj <- tryCatch(
      fit_chains(dat$y, dat$Z, dat$X, iter, warmup, varsel, chain_seeds),
      error = identity
    )
  })[["elapsed"]]
  if (inherits(fit_obj, "error")) {
    row <- error_fit_row(scenario, elapsed, conditionMessage(fit_obj))
  } else {
    row <- summarize_fit(fit_obj, scenario = scenario)
    row$status <- "ok"
    row$error_message <- NA_character_
    write_csv(
      fit_obj$diagnostics$table,
      file.path(out_dir, "parameter_diagnostics",
                sprintf("%s__%s.csv", scenario, row_id))
    )
  }
  row$row_id <- row_id
  row$iter <- iter
  row$warmup <- warmup
  row$varsel <- varsel
  row$chain_seeds <- paste(chain_seeds, collapse = ";")
  row$label_sum <- dat$label_sum
  row
}

save_env_report <- function(path) {
  d <- packageDescription("bkmr")
  news_path <- system.file("NEWS.md", package = "bkmr")
  sink(path)
  on.exit(sink(), add = TRUE)
  cat("Public bkmr verification environment\n")
  cat("Package path:", find.package("bkmr"), "\n")
  cat("Version:", as.character(packageVersion("bkmr")), "\n")
  cat("Date:", if (is.null(d$Date)) NA else d$Date, "\n")
  cat("RemoteSha:", if (is.null(d$RemoteSha)) NA else d$RemoteSha, "\n")
  cat("Built:", if (is.null(d$Built)) NA else d$Built, "\n\n")
  cat("posterior package path:", find.package("posterior"), "\n")
  cat("posterior version:", as.character(packageVersion("posterior")), "\n\n")
  cat("NEWS path:", if (nzchar(news_path)) news_path else NA, "\n\n")
  if (nzchar(news_path)) {
    cat(paste(readLines(news_path, warn = FALSE), collapse = "\n"))
  }
  cat("\n\nsessionInfo()\n")
  print(sessionInfo())
}

dataset_seed <- 303L
iter_main <- 10000L
warmup_main <- 5000L
iter_long <- 100000L
warmup_long <- 50000L

dat40 <- make_bkmr_dataset(dataset_seed, 40L, 3L, 2)
main_chain_seeds <- sample.int(.Machine$integer.max, 4L)
dat80 <- make_bkmr_dataset(dataset_seed, 80L, 3L, 2)

set.seed(20260420L)
seed_sets_3 <- replicate(3L, sample.int(.Machine$integer.max, 4L))
seed_sets_10 <- replicate(10L, sample.int(.Machine$integer.max, 4L))

save_env_report(file.path(out_dir, "verification_env.txt"))

main_rows <- list()
main_rows[[length(main_rows) + 1L]] <- run_existing_data(
  dat40, "step1_main_exact_varsel_true", "step1_main",
  iter_main, warmup_main, TRUE, main_chain_seeds
)
for (i in seq_len(ncol(seed_sets_3))) {
  main_rows[[length(main_rows) + 1L]] <- run_existing_data(
    dat40, sprintf("step3_seed_%d_exact_varsel_true", i), sprintf("step3_seed_%d", i),
    iter_main, warmup_main, TRUE, seed_sets_3[, i]
  )
}
main_rows[[length(main_rows) + 1L]] <- run_existing_data(
  dat40, "step4_long_exact_varsel_false", "step4_long",
  iter_long, warmup_long, FALSE, main_chain_seeds
)
main_verification <- do.call(rbind, main_rows)
main_verification$dataset_seed <- dataset_seed
main_verification$public_bkmr_version <- as.character(packageVersion("bkmr"))
main_verification$public_bkmr_sha <- packageDescription("bkmr")$RemoteSha
write_csv(main_verification, file.path(out_dir, "main_verification_summary.csv"))

crash_result <- tryCatch({
  set.seed(main_chain_seeds[1])
  bkmr::kmbayes(y = dat80$y, Z = dat80$Z, X = dat80$X, iter = iter_main,
                family = "binomial", varsel = FALSE, verbose = FALSE)
  list(status = "ran_to_completion", message = NA_character_)
}, error = function(e) {
  list(status = "error", message = conditionMessage(e))
})
writeLines(c(
  "Step 2: n = 80 numerical crash check",
  paste("dataset_seed:", dataset_seed),
  paste("chain_seed:", main_chain_seeds[1]),
  paste("iter:", iter_main),
  paste("warmup:", warmup_main),
  paste("status:", crash_result$status),
  paste("message:", crash_result$message)
), file.path(out_dir, "crash_verification.txt"))

seed10_rows <- list()
for (i in seq_len(ncol(seed_sets_10))) {
  seed10_rows[[length(seed10_rows) + 1L]] <- run_existing_data(
    dat40, "public_bkmr_seed10_varsel_true", sprintf("seed10_%02d", i),
    iter_main, warmup_main, TRUE, seed_sets_10[, i]
  )
  seed10_rows[[length(seed10_rows)]]$seed_id <- i
}
seed10_details <- do.call(rbind, seed10_rows)
seed10_ok <- seed10_details$status == "ok"
seed10_completed <- seed10_details[seed10_ok, , drop = FALSE]
seed10_summary <- data.frame(
  scenario = "public_bkmr_seed10_varsel_true",
  model = "exact_public",
  dataset_seed = dataset_seed,
  seed_sweep_seed = 20260420L,
  n_seed_sets = 10L,
  n_completed = sum(seed10_ok),
  crash_rate = mean(!seed10_ok),
  hard_failure_rate = if (sum(seed10_ok) == 0L) NA_real_ else mean(!seed10_completed$converged_hard),
  tuning_issue_rate = if (sum(seed10_ok) == 0L) NA_real_ else mean(seed10_completed$tuning_issue),
  median_max_rhat = if (sum(seed10_ok) == 0L) NA_real_ else stats::median(seed10_completed$max_rhat, na.rm = TRUE),
  worst_max_rhat = if (sum(seed10_ok) == 0L) NA_real_ else max(seed10_completed$max_rhat, na.rm = TRUE),
  median_min_ess = if (sum(seed10_ok) == 0L) NA_real_ else stats::median(seed10_completed$min_ess, na.rm = TRUE),
  worst_min_ess = if (sum(seed10_ok) == 0L) NA_real_ else min(seed10_completed$min_ess, na.rm = TRUE),
  median_min_bulk_ess = if (sum(seed10_ok) == 0L) NA_real_ else stats::median(seed10_completed$min_bulk_ess, na.rm = TRUE),
  worst_min_bulk_ess = if (sum(seed10_ok) == 0L) NA_real_ else min(seed10_completed$min_bulk_ess, na.rm = TRUE),
  median_min_tail_ess = if (sum(seed10_ok) == 0L) NA_real_ else stats::median(seed10_completed$min_tail_ess, na.rm = TRUE),
  worst_min_tail_ess = if (sum(seed10_ok) == 0L) NA_real_ else min(seed10_completed$min_tail_ess, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_csv(seed10_details, file.path(out_dir, "seed10_run_details.csv"))
write_csv(seed10_summary, file.path(out_dir, "seed10_summary.csv"))

dataset_sweep_rows <- list()
for (seed in dataset_seeds(20260514L, 50L)) {
  dataset_sweep_rows[[length(dataset_sweep_rows) + 1L]] <- fit_task(
    "public_bkmr_dataset_sweep_n40_varsel_true", seed, 40L, 3L, 2,
    iter_main, warmup_main, TRUE, fixed_chain_seeds(),
    parameter_dir = file.path(out_dir, "parameter_diagnostics")
  )
  dataset_sweep_rows[[length(dataset_sweep_rows) + 1L]] <- fit_task(
    "public_bkmr_dataset_sweep_n80_varsel_false", seed, 80L, 3L, 2,
    iter_main, warmup_main, FALSE, fixed_chain_seeds(),
    parameter_dir = file.path(out_dir, "parameter_diagnostics")
  )
}
dataset_sweep <- do.call(rbind, dataset_sweep_rows)
dataset_sweep$dataset_sweep_seed <- 20260514L
dataset_summary <- summarize_sweep(dataset_sweep, c("scenario", "dataset_sweep_seed"))
write_csv(dataset_sweep, file.path(out_dir, "dataset_sweep_run_details.csv"))
write_csv(dataset_summary, file.path(out_dir, "dataset_sweep_summary.csv"))

n_sweep_rows <- lapply(c(40L, 80L, 160L, 320L), function(n) {
  fit_task("public_bkmr_n_sweep_varsel_false", dataset_seed, n, 3L, 2,
           iter_main, warmup_main, FALSE, fixed_chain_seeds(),
           parameter_dir = file.path(out_dir, "parameter_diagnostics"))
})
n_sweep <- do.call(rbind, n_sweep_rows)
write_csv(n_sweep, file.path(out_dir, "n_sweep_summary.csv"))

long_rows <- lapply(c(50000L, 100000L), function(iter) {
  run_existing_data(dat40, "public_bkmr_longrun_varsel_false",
                    sprintf("long_%d", iter), iter, iter / 2L, FALSE, fixed_chain_seeds())
})
write_csv(do.call(rbind, long_rows), file.path(out_dir, "longrun_50k_100k_summary.csv"))

long3_rows <- list()
for (i in seq_len(ncol(seed_sets_3))) {
  for (iter in c(50000L, 100000L)) {
    long3_rows[[length(long3_rows) + 1L]] <- run_existing_data(
      dat40, "public_bkmr_longrun_seed3_varsel_false",
      sprintf("seed_%d_iter_%d", i, iter), iter, iter / 2L, FALSE, seed_sets_3[, i]
    )
    long3_rows[[length(long3_rows)]]$seed_id <- i
  }
}
long3 <- do.call(rbind, long3_rows)
long3_summary <- summarize_sweep(long3, c("scenario", "iter", "warmup"))
write_csv(long3, file.path(out_dir, "longrun_seed3_run_details.csv"))
write_csv(long3_summary, file.path(out_dir, "longrun_seed3_summary.csv"))
