benchmark_script <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/freeze_tasks_1_3_benchmark.R"
)
source(benchmark_script)

pkg_root <- load_abkmr_benchmark_context()
output_dir <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/results/2026-05-23-binomial-balanced-n40-longrun",
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_seed <- 303L
chain_seeds <- c(142709412L, 1107285817L, 1907135512L, 672667230L)
iter_grid <- c(50000L, 100000L)

generate_dataset <- function(seed) {
  set.seed(seed)
  dat <- bkmr::SimData(
    n = 40,
    M = 4,
    family = "binomial",
    hfun = 2L,
    beta.true = 0.5,
    ind = 1:2
  )
  list(
    y = dat$y,
    Z = as.matrix(dat$Z),
    X = as.matrix(dat$X),
    label_sum = sum(dat$y)
  )
}

dat <- generate_dataset(dataset_seed)
rows <- list()

for (iter in iter_grid) {
  warmup <- as.integer(iter / 2L)
  fit <- tryCatch(
    fit_exact_chains(
      chains = length(chain_seeds),
      seeds = chain_seeds,
      warmup = warmup,
      y = dat$y,
      Z = dat$Z,
      X = dat$X,
      iter = iter,
      verbose = FALSE,
      varsel = FALSE,
      est.h = FALSE,
      family = "binomial"
    ),
    error = identity
  )

  if (inherits(fit, "error")) {
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_seed = dataset_seed,
      n = 40L,
      hfun = 2L,
      beta_true = 0.5,
      iter = iter,
      warmup = warmup,
      label_sum = dat$label_sum,
      status = "error",
      error_message = conditionMessage(fit),
      elapsed_sec = NA_real_,
      max_rhat = NA_real_,
      min_ess = NA_real_,
      converged_hard = FALSE,
      tuning_issue = NA,
      stringsAsFactors = FALSE
    )
  } else {
    run_sum <- summarize_benchmark_run(
      model = "exact_public",
      nd = NA_real_,
      fit_obj = fit,
      exact = TRUE,
      prespecified_univariate = "z1",
      scenario = sprintf("balanced_bkmr_n40_iter_%s", iter)
    )[, c("elapsed_sec", "max_rhat", "min_ess", "converged_hard", "tuning_issue")]

    rows[[length(rows) + 1L]] <- cbind(
      data.frame(
        dataset_seed = dataset_seed,
        n = 40L,
        hfun = 2L,
        beta_true = 0.5,
        iter = iter,
        warmup = warmup,
        label_sum = dat$label_sum,
        status = "ok",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      ),
      run_sum
    )
  }
}

summary_tbl <- do.call(rbind, rows)
meta <- data.frame(
  generator = "bkmr",
  dataset_seed = dataset_seed,
  chain_seed_1 = chain_seeds[1],
  chain_seed_2 = chain_seeds[2],
  chain_seed_3 = chain_seeds[3],
  chain_seed_4 = chain_seeds[4],
  stringsAsFactors = FALSE
)

utils::write.csv(summary_tbl, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(output_dir, "meta.csv"), row.names = FALSE)

print(summary_tbl)
