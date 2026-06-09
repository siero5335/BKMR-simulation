benchmark_script <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/freeze_tasks_1_3_benchmark.R"
)
source(benchmark_script)

pkg_root <- load_abkmr_benchmark_context()
output_dir <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/results/2026-05-22-binomial-default-bkmr-support-sweep",
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_sweep_seed <- 20260522L
K <- 20L
iter <- 10000L
warmup <- 5000L
chain_seeds <- c(142709412L, 1107285817L, 1907135512L, 672667230L)

set.seed(dataset_sweep_seed)
dataset_seeds <- sample.int(.Machine$integer.max, K)

scenario_grid <- data.frame(
  scenario = c(
    "default_bkmr_n40_varsel_true",
    "default_bkmr_n80_varsel_false"
  ),
  n = c(40L, 80L),
  varsel = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)

generate_dataset <- function(seed, n) {
  set.seed(seed)
  dat <- bkmr::SimData(
    n = n,
    M = 4,
    family = "binomial",
    hfun = 3L,
    beta.true = 2,
    ind = 1:2
  )
  list(
    y = dat$y,
    Z = as.matrix(dat$Z),
    X = as.matrix(dat$X),
    label_sum = sum(dat$y)
  )
}

run_rows <- list()
seed_rows <- list()

for (i in seq_len(nrow(scenario_grid))) {
  sc <- scenario_grid[i, ]
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

    dat <- generate_dataset(ds, sc$n)
    fit_obj <- NULL
    elapsed <- system.time({
      fit_obj <- tryCatch(
        fit_exact_chains(
          chains = length(chain_seeds),
          seeds = chain_seeds,
          warmup = warmup,
          y = dat$y,
          Z = dat$Z,
          X = dat$X,
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
        scenario = sc$scenario,
        dataset_seed = ds,
        n = sc$n,
        varsel = sc$varsel,
        label_sum = dat$label_sum,
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
        scenario = sc$scenario,
        dataset_seed = ds,
        n = sc$n,
        varsel = sc$varsel,
        label_sum = dat$label_sum,
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
      n = df$n[1],
      varsel = df$varsel[1],
      dataset_sweep_seed = dataset_sweep_seed,
      K = nrow(df),
      n_completed = sum(ok),
      crash_rate = mean(df$status == "error"),
      hard_failure_rate = if (sum(ok) == 0) NA_real_ else mean(!comp$converged_hard),
      tuning_issue_rate = if (sum(ok) == 0) NA_real_ else mean(comp$tuning_issue, na.rm = TRUE),
      mean_label_sum = mean(df$label_sum, na.rm = TRUE),
      min_label_sum = min(df$label_sum, na.rm = TRUE),
      max_label_sum = max(df$label_sum, na.rm = TRUE),
      median_max_rhat = if (sum(ok) == 0) NA_real_ else stats::median(comp$max_rhat, na.rm = TRUE),
      worst_max_rhat = if (sum(ok) == 0) NA_real_ else max(comp$max_rhat, na.rm = TRUE),
      median_min_ess = if (sum(ok) == 0) NA_real_ else stats::median(comp$min_ess, na.rm = TRUE),
      worst_min_ess = if (sum(ok) == 0) NA_real_ else min(comp$min_ess, na.rm = TRUE),
      n_hard_pass = if (sum(ok) == 0) 0L else sum(comp$converged_hard),
      stringsAsFactors = FALSE
    )
  })
)

meta <- data.frame(
  generator = "bkmr",
  hfun = 3L,
  beta_true = 2,
  dataset_sweep_seed = dataset_sweep_seed,
  K = K,
  iter = iter,
  warmup = warmup,
  chain_seed_1 = chain_seeds[1],
  chain_seed_2 = chain_seeds[2],
  chain_seed_3 = chain_seeds[3],
  chain_seed_4 = chain_seeds[4],
  stringsAsFactors = FALSE
)

utils::write.csv(run_details, file.path(output_dir, "run_details.csv"), row.names = FALSE)
utils::write.csv(summary_tbl, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(seed_table, file.path(output_dir, "seed_table.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(output_dir, "meta.csv"), row.names = FALSE)

print(summary_tbl)
