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

dataset_sweep_seed <- 20260514L
K <- 50L
iter <- 10000L
warmup <- 5000L
chain_seeds <- c(142709412L, 1107285817L, 1907135512L, 672667230L)

set.seed(dataset_sweep_seed)
dataset_seeds <- sample.int(.Machine$integer.max, K)

scenarios <- list(
  list(
    scenario = "public_bkmr_dataset_sweep_n40_varsel_true",
    n = 40L,
    varsel = TRUE
  ),
  list(
    scenario = "public_bkmr_dataset_sweep_n80_varsel_false",
    n = 80L,
    varsel = FALSE
  )
)

run_rows <- list()
seed_rows <- list()

for (sc in scenarios) {
  for (ds in dataset_seeds) {
    seed_rows[[length(seed_rows) + 1L]] <- data.frame(
      scenario = sc$scenario,
      dataset_seed = ds,
      chain_1 = chain_seeds[1],
      chain_2 = chain_seeds[2],
      chain_3 = chain_seeds[3],
      chain_4 = chain_seeds[4],
      stringsAsFactors = FALSE
    )

    set.seed(ds)
    dat <- sim_data(n = sc$n, M = 4, family = "binomial", ind = 1:2)
    Z <- as.matrix(dat$Z)
    X <- as.matrix(dat$X)
    y <- dat$y

    fit_obj <- NULL
    elapsed <- system.time({
      fit_obj <- tryCatch(
        fit_exact_chains(
          chains = length(chain_seeds),
          seeds = chain_seeds,
          warmup = warmup,
          y = y,
          Z = Z,
          X = X,
          iter = iter,
          family = "binomial",
          varsel = sc$varsel,
          verbose = FALSE,
          est.h = FALSE
        ),
        error = identity
      )
    })[["elapsed"]]

    if (inherits(fit_obj, "error")) {
      run_rows[[length(run_rows) + 1L]] <- data.frame(
        dataset_seed = ds,
        scenario = sc$scenario,
        n = sc$n,
        varsel = sc$varsel,
        iter = iter,
        warmup = warmup,
        status = "error",
        error_message = conditionMessage(fit_obj),
        max_rhat = NA_real_,
        min_ess = NA_real_,
        converged_hard = FALSE,
        tuning_issue = NA,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
    } else {
      row <- summarize_benchmark_run(
        model = "exact_public",
        nd = NA_real_,
        fit_obj = fit_obj,
        exact = TRUE,
        prespecified_univariate = "z1",
        scenario = sc$scenario
      )
      run_rows[[length(run_rows) + 1L]] <- data.frame(
        dataset_seed = ds,
        scenario = sc$scenario,
        n = sc$n,
        varsel = sc$varsel,
        iter = iter,
        warmup = warmup,
        status = "ok",
        error_message = NA_character_,
        max_rhat = row$max_rhat,
        min_ess = row$min_ess,
        converged_hard = row$converged_hard,
        tuning_issue = row$tuning_issue,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
    }
  }
}

run_details <- do.call(rbind, run_rows)
seed_table <- do.call(rbind, seed_rows)

summary_tbl <- do.call(
  rbind,
  lapply(split(run_details, run_details$scenario), function(df) {
    ok <- df$status == "ok"
    comp <- df[ok, , drop = FALSE]
    data.frame(
      scenario = df$scenario[1],
      dataset_sweep_seed = dataset_sweep_seed,
      K = nrow(df),
      n_completed = sum(ok),
      crash_rate = mean(df$status == "error"),
      hard_failure_rate = if (sum(ok) == 0) NA_real_ else mean(!comp$converged_hard),
      tuning_issue_rate = if (sum(ok) == 0) NA_real_ else mean(comp$tuning_issue, na.rm = TRUE),
      median_max_rhat = if (sum(ok) == 0) NA_real_ else stats::median(comp$max_rhat, na.rm = TRUE),
      worst_max_rhat = if (sum(ok) == 0) NA_real_ else max(comp$max_rhat, na.rm = TRUE),
      median_min_ess = if (sum(ok) == 0) NA_real_ else stats::median(comp$min_ess, na.rm = TRUE),
      worst_min_ess = if (sum(ok) == 0) NA_real_ else min(comp$min_ess, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

utils::write.csv(
  run_details,
  file.path(results_dir, "dataset_sweep_run_details.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_tbl,
  file.path(results_dir, "dataset_sweep_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  seed_table,
  file.path(results_dir, "dataset_sweep_seed_table.csv"),
  row.names = FALSE
)

print(summary_tbl)
