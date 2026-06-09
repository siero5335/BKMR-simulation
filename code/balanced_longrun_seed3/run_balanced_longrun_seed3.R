benchmark_script <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/freeze_tasks_1_3_benchmark.R"
)
source(benchmark_script)

pkg_root <- load_abkmr_benchmark_context()
output_dir <- normalizePath(
  "/Users/a_eguchi/Documents/GitHub/abkmr-dev/A-BKMR-main/inst/benchmarks/results/2026-05-23-binomial-balanced-longrun-seed3",
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_sweep_seed <- 20260522L
chain_seeds <- c(142709412L, 1107285817L, 1907135512L, 672667230L)
iter_grid <- c(50000L, 100000L)

set.seed(dataset_sweep_seed)
dataset_seeds <- sample.int(.Machine$integer.max, 3L)

scenarios <- list(
  list(
    scenario = "balanced_bkmr_n40_varsel_false_longrun",
    n = 40L,
    varsel = FALSE
  ),
  list(
    scenario = "balanced_bkmr_n80_varsel_false_longrun",
    n = 80L,
    varsel = FALSE
  )
)

generate_dataset <- function(seed, n) {
  set.seed(seed)
  dat <- bkmr::SimData(
    n = n,
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

rows <- list()
seed_rows <- list()

for (sc in scenarios) {
  for (dataset_seed in dataset_seeds) {
    dat <- generate_dataset(dataset_seed, sc$n)

    for (iter in iter_grid) {
      warmup <- as.integer(iter / 2L)

      seed_rows[[length(seed_rows) + 1L]] <- data.frame(
        scenario = sc$scenario,
        n = sc$n,
        varsel = sc$varsel,
        dataset_seed = dataset_seed,
        iter = iter,
        warmup = warmup,
        chain_1 = chain_seeds[1],
        chain_2 = chain_seeds[2],
        chain_3 = chain_seeds[3],
        chain_4 = chain_seeds[4],
        stringsAsFactors = FALSE
      )

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
          varsel = sc$varsel,
          est.h = FALSE,
          family = "binomial"
        ),
        error = identity
      )

      if (inherits(fit, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          scenario = sc$scenario,
          dataset_seed = dataset_seed,
          n = sc$n,
          varsel = sc$varsel,
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
          scenario = sprintf("%s_iter_%s", sc$scenario, iter)
        )[, c("elapsed_sec", "max_rhat", "min_ess", "converged_hard", "tuning_issue")]

        rows[[length(rows) + 1L]] <- cbind(
          data.frame(
            scenario = sc$scenario,
            dataset_seed = dataset_seed,
            n = sc$n,
            varsel = sc$varsel,
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
  }
}

run_details <- do.call(rbind, rows)
seed_table <- do.call(rbind, seed_rows)

summary_tbl <- do.call(
  rbind,
  lapply(split(run_details, list(run_details$scenario, run_details$iter), drop = TRUE), function(df) {
    ok <- df$status == "ok"
    comp <- df[ok, , drop = FALSE]
    data.frame(
      scenario = df$scenario[1],
      n = df$n[1],
      varsel = df$varsel[1],
      hfun = df$hfun[1],
      beta_true = df$beta_true[1],
      iter = df$iter[1],
      warmup = df$warmup[1],
      K = nrow(df),
      n_completed = sum(ok),
      crash_rate = mean(df$status == "error"),
      hard_failure_rate = if (sum(ok) == 0) NA_real_ else mean(!comp$converged_hard),
      tuning_issue_rate = if (sum(ok) == 0) NA_real_ else mean(comp$tuning_issue, na.rm = TRUE),
      median_label_sum = stats::median(df$label_sum),
      median_max_rhat = if (sum(ok) == 0) NA_real_ else stats::median(comp$max_rhat, na.rm = TRUE),
      worst_max_rhat = if (sum(ok) == 0) NA_real_ else max(comp$max_rhat, na.rm = TRUE),
      median_min_ess = if (sum(ok) == 0) NA_real_ else stats::median(comp$min_ess, na.rm = TRUE),
      worst_min_ess = if (sum(ok) == 0) NA_real_ else min(comp$min_ess, na.rm = TRUE),
      n_hard_pass = if (sum(ok) == 0) 0L else sum(comp$converged_hard, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

meta <- data.frame(
  generator = "bkmr",
  dataset_sweep_seed = dataset_sweep_seed,
  dataset_seed_1 = dataset_seeds[1],
  dataset_seed_2 = dataset_seeds[2],
  dataset_seed_3 = dataset_seeds[3],
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
