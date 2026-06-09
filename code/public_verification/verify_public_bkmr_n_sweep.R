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
chain_seeds <- c(142709412L, 1107285817L, 1907135512L, 672667230L)
iter <- 10000L
warmup <- 5000L
n_grid <- c(40L, 80L, 160L, 320L)

run_rows <- list()

for (n_val in n_grid) {
  set.seed(dataset_seed)
  dat <- sim_data(n = n_val, M = 4, family = "binomial", ind = 1:2)
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
        varsel = FALSE,
        verbose = FALSE,
        est.h = FALSE
      ),
      error = identity
    )
  })[["elapsed"]]

  if (inherits(fit_obj, "error")) {
    run_rows[[length(run_rows) + 1L]] <- data.frame(
      scenario = "public_bkmr_n_sweep_varsel_false",
      model = "exact_public",
      n = n_val,
      iter = iter,
      warmup = warmup,
      elapsed_sec = elapsed,
      status = "error",
      error_message = conditionMessage(fit_obj),
      max_rhat = NA_real_,
      min_ess = NA_real_,
      converged_hard = FALSE,
      tuning_issue = NA,
      stringsAsFactors = FALSE
    )
  } else {
    row <- summarize_benchmark_run(
      model = "exact_public",
      nd = NA_real_,
      fit_obj = fit_obj,
      exact = TRUE,
      prespecified_univariate = "z1",
      scenario = "public_bkmr_n_sweep_varsel_false"
    )
    row$n <- n_val
    row$iter <- iter
    row$warmup <- warmup
    row$status <- "ok"
    row$error_message <- NA_character_
    run_rows[[length(run_rows) + 1L]] <- row[, c(
      "scenario", "model", "n", "iter", "warmup", "elapsed_sec", "status",
      "error_message", "max_rhat", "min_ess", "converged_hard", "tuning_issue"
    )]
  }
}

n_sweep <- do.call(rbind, run_rows)
meta <- data.frame(
  dataset_seed = dataset_seed,
  chain_seeds = paste(chain_seeds, collapse = ";"),
  iter = iter,
  warmup = warmup,
  public_bkmr_version = as.character(packageVersion("bkmr")),
  public_bkmr_sha = packageDescription("bkmr")$RemoteSha,
  stringsAsFactors = FALSE
)

utils::write.csv(n_sweep, file.path(results_dir, "n_sweep_summary.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(results_dir, "n_sweep_meta.csv"), row.names = FALSE)

print(n_sweep)
