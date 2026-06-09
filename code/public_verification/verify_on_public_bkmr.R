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
  cat("NEWS path:", if (nzchar(news_path)) news_path else NA, "\n\n")
  cat("NEWS excerpt\n")
  cat("===========\n")
  if (nzchar(news_path)) {
    cat(paste(readLines(news_path, warn = FALSE), collapse = "\n"))
  }
  cat("\n\nsessionInfo()\n")
  cat("=============\n")
  print(sessionInfo())
}

capture_fit_summary <- function(model, fit_obj, exact, scenario, nd,
                                prespecified_univariate = "z1") {
  row <- summarize_benchmark_run(
    model = model,
    nd = nd,
    fit_obj = fit_obj,
    exact = exact,
    prespecified_univariate = prespecified_univariate,
    scenario = scenario
  )
  row$status <- "ok"
  row$error_message <- NA_character_
  row
}

capture_fit_error <- function(model, scenario, nd, elapsed, error_message) {
  data.frame(
    scenario = scenario,
    model = model,
    nd = nd,
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
    error_message = as.character(error_message),
    stringsAsFactors = FALSE
  )
}

run_exact_multichain <- function(y, Z, X, iter, warmup, varsel, chain_seeds,
                                 scenario, row_id) {
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
        varsel = varsel,
        verbose = FALSE,
        est.h = FALSE
      ),
      error = identity
    )
  })[["elapsed"]]

  out <- if (inherits(fit_obj, "error")) {
    capture_fit_error(
      model = "exact_public",
      scenario = scenario,
      nd = NA_real_,
      elapsed = elapsed,
      error_message = conditionMessage(fit_obj)
    )
  } else {
    row <- capture_fit_summary(
      model = "exact_public",
      fit_obj = fit_obj,
      exact = TRUE,
      scenario = scenario,
      nd = NA_real_
    )
    row$elapsed_sec <- elapsed
    row
  }

  out$row_id <- row_id
  out$iter <- iter
  out$warmup <- warmup
  out$varsel <- varsel
  out$chain_seeds <- paste(chain_seeds, collapse = ";")
  out
}

dataset_seed <- 303L
chains <- 4L
iter_main <- 10000L
warmup_main <- 5000L
iter_long <- 100000L
warmup_long <- 50000L

set.seed(dataset_seed)
dat40 <- sim_data(n = 40, M = 4, family = "binomial", ind = 1:2)
Z40 <- as.matrix(dat40$Z)
X40 <- as.matrix(dat40$X)
y40 <- dat40$y
main_chain_seeds <- sample.int(.Machine$integer.max, chains)

set.seed(dataset_seed)
dat80 <- sim_data(n = 80, M = 4, family = "binomial", ind = 1:2)
Z80 <- as.matrix(dat80$Z)
X80 <- as.matrix(dat80$X)
y80 <- dat80$y

seed_robustness_seed <- 20260420L
set.seed(seed_robustness_seed)
robustness_seeds <- replicate(3, sample.int(.Machine$integer.max, chains))

save_env_report(file.path(results_dir, "verification_env.txt"))

main_rows <- list()
seed_rows <- list()

main_rows[[1]] <- run_exact_multichain(
  y = y40,
  Z = Z40,
  X = X40,
  iter = iter_main,
  warmup = warmup_main,
  varsel = TRUE,
  chain_seeds = main_chain_seeds,
  scenario = "step1_main_exact_varsel_true",
  row_id = "step1_main"
)

for (i in seq_len(ncol(robustness_seeds))) {
  main_rows[[length(main_rows) + 1L]] <- run_exact_multichain(
    y = y40,
    Z = Z40,
    X = X40,
    iter = iter_main,
    warmup = warmup_main,
    varsel = TRUE,
    chain_seeds = robustness_seeds[, i],
    scenario = sprintf("step3_seed_%d_exact_varsel_true", i),
    row_id = sprintf("step3_seed_%d", i)
  )
}

main_rows[[length(main_rows) + 1L]] <- run_exact_multichain(
  y = y40,
  Z = Z40,
  X = X40,
  iter = iter_long,
  warmup = warmup_long,
  varsel = FALSE,
  chain_seeds = main_chain_seeds,
  scenario = "step4_long_exact_varsel_false",
  row_id = "step4_long"
)

main_verification <- do.call(rbind, main_rows)
main_verification$dataset_seed <- dataset_seed
main_verification$public_bkmr_version <- as.character(packageVersion("bkmr"))
main_verification$public_bkmr_sha <- packageDescription("bkmr")$RemoteSha

crash_result <- tryCatch({
  set.seed(main_chain_seeds[1])
  bkmr::kmbayes(
    y = y80,
    Z = Z80,
    X = X80,
    iter = iter_main,
    family = "binomial",
    varsel = FALSE,
    verbose = FALSE
  )
  list(status = "ran_to_completion", message = NA_character_)
}, error = function(e) {
  list(status = "error", message = conditionMessage(e))
})

crash_lines <- c(
  "Step 2: n = 80 numerical crash check",
  paste("dataset_seed:", dataset_seed),
  paste("chain_seed:", main_chain_seeds[1]),
  paste("iter:", iter_main),
  paste("warmup:", warmup_main),
  paste("status:", crash_result$status),
  paste("message:", if (is.null(crash_result$message)) NA_character_ else crash_result$message)
)
writeLines(crash_lines, file.path(results_dir, "crash_verification.txt"))

seed_table <- do.call(rbind, c(
  list(data.frame(
    label = "step1_main",
    chain_1 = main_chain_seeds[1],
    chain_2 = main_chain_seeds[2],
    chain_3 = main_chain_seeds[3],
    chain_4 = main_chain_seeds[4],
    stringsAsFactors = FALSE
  )),
  lapply(seq_len(ncol(robustness_seeds)), function(i) {
    data.frame(
      label = sprintf("step3_seed_%d", i),
      chain_1 = robustness_seeds[1, i],
      chain_2 = robustness_seeds[2, i],
      chain_3 = robustness_seeds[3, i],
      chain_4 = robustness_seeds[4, i],
      stringsAsFactors = FALSE
    )
  }),
  list(data.frame(
    label = "step4_long",
    chain_1 = main_chain_seeds[1],
    chain_2 = main_chain_seeds[2],
    chain_3 = main_chain_seeds[3],
    chain_4 = main_chain_seeds[4],
    stringsAsFactors = FALSE
  ))
))

meta <- data.frame(
  dataset_seed = dataset_seed,
  seed_robustness_seed = seed_robustness_seed,
  n_main = nrow(Z40),
  n_crash = nrow(Z80),
  iter_main = iter_main,
  warmup_main = warmup_main,
  iter_long = iter_long,
  warmup_long = warmup_long,
  public_bkmr_version = as.character(packageVersion("bkmr")),
  public_bkmr_sha = packageDescription("bkmr")$RemoteSha,
  public_bkmr_path = find.package("bkmr"),
  stringsAsFactors = FALSE
)

utils::write.csv(main_verification, file.path(results_dir, "main_verification_summary.csv"), row.names = FALSE)
utils::write.csv(seed_table, file.path(results_dir, "seed_table.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(results_dir, "meta.csv"), row.names = FALSE)

print(main_verification[, c(
  "row_id", "iter", "warmup", "varsel", "status",
  "max_rhat", "min_ess", "converged_hard", "tuning_issue"
)])
