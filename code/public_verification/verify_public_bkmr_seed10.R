pkg_root <- "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main"
results_dir <- file.path(
  pkg_root,
  "inst/benchmarks/results/2026-04-20-public-bkmr-verification-bundle"
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(RETICULATE_PYTHON = file.path(pkg_root, ".venv/bin/python"))
.libPaths(c(.libPaths(), file.path(pkg_root, ".r-lib")))

library(bkmr, lib.loc = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library")
pkgload::load_all(pkg_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)
source(file.path(pkg_root, "inst/benchmarks/freeze_tasks_1_3_benchmark.R"))

dataset_seed <- 303L
seed_sweep_seed <- 20260420L
n_seed_sets <- 10L
chains <- 4L
iter <- 10000L
warmup <- 5000L

set.seed(dataset_seed)
dat40 <- sim_data(n = 40, M = 4, family = "binomial", ind = 1:2)
Z40 <- as.matrix(dat40$Z)
X40 <- as.matrix(dat40$X)
y40 <- dat40$y

set.seed(seed_sweep_seed)
seed_sets <- replicate(n_seed_sets, sample.int(.Machine$integer.max, chains))

run_rows <- list()
seed_rows <- list()

for (i in seq_len(ncol(seed_sets))) {
  chain_seeds <- seed_sets[, i]
  seed_rows[[length(seed_rows) + 1L]] <- data.frame(
    label = sprintf("seed10_%02d", i),
    chain_1 = chain_seeds[1],
    chain_2 = chain_seeds[2],
    chain_3 = chain_seeds[3],
    chain_4 = chain_seeds[4],
    stringsAsFactors = FALSE
  )

  fit_obj <- NULL
  elapsed <- system.time({
    fit_obj <- tryCatch(
      fit_exact_chains(
        chains = chains,
        seeds = chain_seeds,
        warmup = warmup,
        y = y40,
        Z = Z40,
        X = X40,
        iter = iter,
        family = "binomial",
        varsel = TRUE,
        verbose = FALSE,
        est.h = FALSE
      ),
      error = identity
    )
  })[["elapsed"]]

  if (inherits(fit_obj, "error")) {
    run_rows[[length(run_rows) + 1L]] <- data.frame(
      scenario = "public_bkmr_seed10_varsel_true",
      model = "exact_public",
      nd = NA_real_,
      elapsed_sec = elapsed,
      median_ess_per_sec = NA_real_,
      converged_hard = FALSE,
      tuning_issue = NA,
      converged_strict = FALSE,
      converged = FALSE,
      max_rhat = NA_real_,
      min_ess = NA_real_,
      max_mcse = NA_real_,
      n_rhat_flag = NA_integer_,
      n_ess_flag = NA_integer_,
      n_acceptance_flag = NA_integer_,
      n_hard_flag = NA_integer_,
      n_flagged = NA_integer_,
      overall_q25_to_q75 = NA_real_,
      prespecified_univar_var = "z1",
      prespecified_univar_est = NA_real_,
      top_pip_var = NA_character_,
      top_pip = NA_real_,
      status = "error",
      error_message = conditionMessage(fit_obj),
      seed_id = i,
      iter = iter,
      warmup = warmup,
      varsel = TRUE,
      chain_seeds = paste(chain_seeds, collapse = ";"),
      stringsAsFactors = FALSE
    )
  } else {
    row <- summarize_benchmark_run(
      model = "exact_public",
      nd = NA_real_,
      fit_obj = fit_obj,
      exact = TRUE,
      prespecified_univariate = "z1",
      scenario = "public_bkmr_seed10_varsel_true"
    )
    row$elapsed_sec <- elapsed
    row$status <- "ok"
    row$error_message <- NA_character_
    row$seed_id <- i
    row$iter <- iter
    row$warmup <- warmup
    row$varsel <- TRUE
    row$chain_seeds <- paste(chain_seeds, collapse = ";")
    run_rows[[length(run_rows) + 1L]] <- row
  }
}

run_details <- do.call(rbind, run_rows)
seed_table <- do.call(rbind, seed_rows)
ok <- run_details$status == "ok"
completed <- run_details[ok, , drop = FALSE]

summary <- data.frame(
  scenario = "public_bkmr_seed10_varsel_true",
  model = "exact_public",
  dataset_seed = dataset_seed,
  seed_sweep_seed = seed_sweep_seed,
  n_seed_sets = n_seed_sets,
  n_completed = sum(ok),
  crash_rate = mean(!ok),
  hard_failure_rate = if (sum(ok) == 0) NA_real_ else mean(!completed$converged_hard),
  tuning_issue_rate = if (sum(ok) == 0) NA_real_ else mean(completed$tuning_issue),
  median_max_rhat = if (sum(ok) == 0) NA_real_ else stats::median(completed$max_rhat, na.rm = TRUE),
  worst_max_rhat = if (sum(ok) == 0) NA_real_ else max(completed$max_rhat, na.rm = TRUE),
  median_min_ess = if (sum(ok) == 0) NA_real_ else stats::median(completed$min_ess, na.rm = TRUE),
  worst_min_ess = if (sum(ok) == 0) NA_real_ else min(completed$min_ess, na.rm = TRUE),
  stringsAsFactors = FALSE
)

utils::write.csv(run_details, file.path(results_dir, "seed10_run_details.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(results_dir, "seed10_summary.csv"), row.names = FALSE)
utils::write.csv(seed_table, file.path(results_dir, "seed10_seed_table.csv"), row.names = FALSE)

print(summary)
