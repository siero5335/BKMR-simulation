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
n_seed_sets <- 3L
chains <- 4L
configs <- list(
  list(iter = 50000L, warmup = 25000L),
  list(iter = 100000L, warmup = 50000L)
)

set.seed(dataset_seed)
dat <- sim_data(n = 40, M = 4, family = "binomial", ind = 1:2)
Z <- as.matrix(dat$Z)
X <- as.matrix(dat$X)
y <- dat$y

set.seed(seed_sweep_seed)
seed_sets <- replicate(n_seed_sets, sample.int(.Machine$integer.max, chains))

run_rows <- list()
seed_rows <- list()

for (seed_id in seq_len(ncol(seed_sets))) {
  chain_seeds <- seed_sets[, seed_id]
  seed_rows[[length(seed_rows) + 1L]] <- data.frame(
    seed_id = seed_id,
    chain_1 = chain_seeds[1],
    chain_2 = chain_seeds[2],
    chain_3 = chain_seeds[3],
    chain_4 = chain_seeds[4],
    stringsAsFactors = FALSE
  )

  for (cfg in configs) {
    fit_obj <- NULL
    elapsed <- system.time({
      fit_obj <- tryCatch(
        fit_exact_chains(
          chains = chains,
          seeds = chain_seeds,
          warmup = cfg$warmup,
          y = y,
          Z = Z,
          X = X,
          iter = cfg$iter,
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
        scenario = "public_bkmr_longrun_seed3_varsel_false",
        model = "exact_public",
        seed_id = seed_id,
        iter = cfg$iter,
        warmup = cfg$warmup,
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
        scenario = "public_bkmr_longrun_seed3_varsel_false"
      )
      row$seed_id <- seed_id
      row$iter <- cfg$iter
      row$warmup <- cfg$warmup
      row$status <- "ok"
      row$error_message <- NA_character_
      run_rows[[length(run_rows) + 1L]] <- row[, c(
        "scenario", "model", "seed_id", "iter", "warmup", "elapsed_sec", "status",
        "error_message", "max_rhat", "min_ess", "converged_hard", "tuning_issue"
      )]
    }
  }
}

run_details <- do.call(rbind, run_rows)
seed_table <- do.call(rbind, seed_rows)

summary <- do.call(
  rbind,
  lapply(split(run_details, run_details$iter), function(df) {
    ok <- df$status == "ok"
    completed <- df[ok, , drop = FALSE]
    data.frame(
      scenario = df$scenario[1],
      model = df$model[1],
      iter = df$iter[1],
      warmup = df$warmup[1],
      n_seed_sets = nrow(df),
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
  })
)

meta <- data.frame(
  dataset_seed = dataset_seed,
  seed_sweep_seed = seed_sweep_seed,
  n_seed_sets = n_seed_sets,
  n = nrow(Z),
  M = ncol(Z),
  varsel = FALSE,
  public_bkmr_version = as.character(packageVersion("bkmr")),
  public_bkmr_sha = packageDescription("bkmr")$RemoteSha,
  stringsAsFactors = FALSE
)

utils::write.csv(run_details, file.path(results_dir, "longrun_seed3_run_details.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(results_dir, "longrun_seed3_summary.csv"), row.names = FALSE)
utils::write.csv(seed_table, file.path(results_dir, "longrun_seed3_seed_table.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(results_dir, "longrun_seed3_meta.csv"), row.names = FALSE)

print(summary)
