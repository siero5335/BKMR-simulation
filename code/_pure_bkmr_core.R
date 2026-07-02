ensure_bkmr_available <- function() {
  if (!requireNamespace("bkmr", quietly = TRUE)) {
    stop("The bkmr package is required.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("The posterior package is required for Vehtari et al. diagnostics.", call. = FALSE)
  }
  invisible(TRUE)
}

project_root <- function(script_file) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
}

script_path <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else ""
  }, character(1))
  if (any(nzchar(ofiles))) {
    return(normalizePath(tail(ofiles[nzchar(ofiles)], 1), mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Could not locate the current script path.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

safe_file_component <- function(x) {
  gsub("[^A-Za-z0-9_.=-]+", "_", as.character(x))
}

make_bkmr_dataset <- function(seed, n, hfun, beta_true, M = 4L, ind = 1:2) {
  set.seed(seed)
  dat <- bkmr::SimData(
    n = n,
    M = M,
    family = "binomial",
    hfun = hfun,
    beta.true = beta_true,
    ind = ind
  )
  list(
    y = dat$y,
    Z = as.matrix(dat$Z),
    X = as.matrix(dat$X),
    label_sum = sum(dat$y)
  )
}

post_warmup_index <- function(fit, warmup = NULL) {
  if (is.null(warmup)) {
    warmup <- floor(fit$iter / 2)
  }
  keep <- seq.int(min(fit$iter, warmup + 1L), fit$iter)
  if (length(keep) < 4L) {
    stop("Need at least four retained iterations per chain.", call. = FALSE)
  }
  keep
}

param_label <- function(prefix, idx, labels = NULL) {
  if (!is.null(labels) && length(labels) >= idx && !is.na(labels[idx]) && nzchar(labels[idx])) {
    return(sprintf("%s[%s]", prefix, labels[idx]))
  }
  if (identical(prefix, "lambda") && idx == 1L) {
    return("lambda")
  }
  sprintf("%s[%d]", prefix, idx)
}

param_draws <- function(chains, extractor, warmup) {
  keep <- lapply(chains, post_warmup_index, warmup = warmup)
  kept <- vapply(keep, length, integer(1))
  if (length(unique(kept)) != 1L) {
    stop("All chains must have the same number of retained iterations.", call. = FALSE)
  }
  vals <- lapply(seq_along(chains), function(i) as.numeric(extractor(chains[[i]])[keep[[i]]]))
  do.call(cbind, vals)
}

constant_draws <- function(draws) {
  vals <- as.vector(draws)
  vals <- vals[is.finite(vals)]
  length(vals) > 0L && isTRUE(all(vals == vals[1]))
}

bkmr_starting_values_description <- function() {
  paste(
    "chain-specific randomized binomial starts from fixed initial-value seed;",
    "h.hat ~ Uniform(-2, 2) and beta ~ Uniform(-1, 1);",
    "positive r and |ystar| generated as exp(Uniform(-1, 1));",
    "lambda generated as log-Uniform(1, 20);",
    "ystar sign constrained by y; sigsq.eps fixed at 1;",
    "delta ~ Bernoulli(0.5) for varsel=TRUE with r=0 when delta=0"
  )
}

fixed_initial_value_seed <- function() {
  20260621L
}

derive_initial_value_seed <- function(base_seed = fixed_initial_value_seed(), ...) {
  modulus <- .Machine$integer.max
  key <- paste(c(as.character(base_seed), as.character(unlist(list(...), use.names = FALSE))),
               collapse = "|")
  seed <- as.numeric(base_seed %% modulus)
  for (code in utf8ToInt(key)) {
    seed <- (seed * 131 + code) %% modulus
  }
  seed <- as.integer(seed)
  if (is.na(seed) || seed <= 0L) {
    seed <- 1L
  }
  seed
}

with_preserved_seed <- function(seed, expr) {
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

bkmr_binomial_starting_values <- function(y, Z, X, varsel) {
  delta <- if (isTRUE(varsel)) stats::rbinom(ncol(Z), size = 1L, prob = 0.5) else rep(1, ncol(Z))
  r <- exp(stats::runif(ncol(Z), min = -1, max = 1))
  if (isTRUE(varsel)) {
    r[delta == 0L] <- 0
  }
  list(
    h.hat = stats::runif(length(y), min = -2, max = 2),
    beta = stats::runif(ncol(X), min = -1, max = 1),
    sigsq.eps = 1,
    r = r,
    lambda = exp(stats::runif(1L, min = log(1), max = log(20))),
    delta = delta,
    ystar = ifelse(y == 1, 1, -1) * exp(stats::runif(length(y), min = -1, max = 1))
  )
}

make_bkmr_binomial_starting_values <- function(y, Z, X, n_chains, varsel,
                                               initial_value_seed = fixed_initial_value_seed()) {
  with_preserved_seed(initial_value_seed, {
    lapply(seq_len(n_chains), function(...) {
      bkmr_binomial_starting_values(y, Z, X, varsel = varsel)
    })
  })
}

posterior_stat <- function(draws, fun) {
  tryCatch(as.numeric(fun(draws)), error = function(...) NA_real_)
}

diagnostic_draw_table <- function(chains, warmup) {
  rows <- list()
  beta_labels <- colnames(chains[[1]]$X)
  z_labels <- colnames(chains[[1]]$Z)

  rows[[length(rows) + 1L]] <- data.frame(
    parameter = "sigsq.eps",
    draws = I(list(param_draws(chains, function(fit) fit$sigsq.eps, warmup)))
  )

  lambda <- as.matrix(chains[[1]]$lambda)
  for (j in seq_len(ncol(lambda))) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = param_label("lambda", j),
      draws = I(list(param_draws(chains, function(fit) as.matrix(fit$lambda)[, j], warmup)))
    )
  }

  r <- as.matrix(chains[[1]]$r)
  for (j in seq_len(ncol(r))) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = param_label("r", j, z_labels),
      draws = I(list(param_draws(chains, function(fit) as.matrix(fit$r)[, j], warmup)))
    )
  }

  beta <- as.matrix(chains[[1]]$beta)
  for (j in seq_len(ncol(beta))) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = param_label("beta", j, beta_labels),
      draws = I(list(param_draws(chains, function(fit) as.matrix(fit$beta)[, j], warmup)))
    )
  }

  do.call(rbind, rows)
}

acceptance_table <- function(chains, warmup) {
  rows <- list()
  lambda_acc <- as.matrix(chains[[1]]$acc.lambda)
  for (j in seq_len(ncol(lambda_acc))) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = param_label("lambda", j),
      acceptance = mean(vapply(chains, function(fit) {
        acc <- as.matrix(fit$acc.lambda)
        mean(acc[post_warmup_index(fit, warmup), j], na.rm = TRUE)
      }, numeric(1))),
      stringsAsFactors = FALSE
    )
  }

  r_acc <- as.matrix(chains[[1]]$acc.r)
  for (j in seq_len(ncol(r_acc))) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = param_label("r", j, colnames(chains[[1]]$Z)),
      acceptance = mean(vapply(chains, function(fit) {
        acc <- as.matrix(fit$acc.r)
        mean(acc[post_warmup_index(fit, warmup), j], na.rm = TRUE)
      }, numeric(1))),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(chains[[1]]$acc.rdelta)) {
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = "rdelta_move",
      acceptance = mean(vapply(chains, function(fit) {
        mean(fit$acc.rdelta[post_warmup_index(fit, warmup)], na.rm = TRUE)
      }, numeric(1))),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

diagnose_chains <- function(chains, warmup, rhat_threshold = 1.01,
                            ess_threshold = 400,
                            acceptance_bounds = c(0.1, 0.8)) {
  draw_rows <- diagnostic_draw_table(chains, warmup)
  table <- do.call(rbind, lapply(seq_len(nrow(draw_rows)), function(i) {
    draws <- draw_rows$draws[[i]]
    is_constant <- constant_draws(draws)
    data.frame(
      parameter = draw_rows$parameter[i],
      rhat = posterior_stat(draws, posterior::rhat),
      bulk_ess = posterior_stat(draws, posterior::ess_bulk),
      tail_ess = posterior_stat(draws, posterior::ess_tail),
      mcse_mean = posterior_stat(draws, posterior::mcse_mean),
      constant_draws = is_constant,
      diagnostic_note = if (is_constant) {
        "constant/fixed draws; posterior R-hat and ESS are undefined"
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }))

  acc <- acceptance_table(chains, warmup)
  table <- merge(table, acc, by = "parameter", all = TRUE, sort = FALSE)
  table$rhat_flag <- !is.na(table$rhat) & table$rhat > rhat_threshold
  table$bulk_ess_flag <- !is.na(table$bulk_ess) & table$bulk_ess < ess_threshold
  table$tail_ess_flag <- !is.na(table$tail_ess) & table$tail_ess < ess_threshold
  table$ess_flag <- table$bulk_ess_flag | table$tail_ess_flag
  table$acceptance_flag <- !is.na(table$acceptance) &
    (table$acceptance < acceptance_bounds[1] | table$acceptance > acceptance_bounds[2])
  table$hard_flag <- table$rhat_flag | table$ess_flag
  table$tuning_flag <- table$acceptance_flag
  table$needs_attention <- table$hard_flag | table$tuning_flag

  list(table = table, acceptance = acc, warmup = warmup)
}

fit_chains <- function(y, Z, X, iter, warmup, varsel, chain_seeds,
                       initial_value_seed = fixed_initial_value_seed()) {
  chains <- vector("list", length(chain_seeds))
  starting_values <- make_bkmr_binomial_starting_values(
    y = y,
    Z = Z,
    X = X,
    n_chains = length(chain_seeds),
    varsel = varsel,
    initial_value_seed = initial_value_seed
  )
  elapsed <- system.time({
    for (i in seq_along(chain_seeds)) {
      set.seed(chain_seeds[i])
      invisible(utils::capture.output({
        chains[[i]] <- suppressMessages(bkmr::kmbayes(
          y = y,
          Z = Z,
          X = X,
          iter = iter,
          family = "binomial",
          varsel = varsel,
          verbose = FALSE,
          est.h = FALSE,
          starting.values = starting_values[[i]]
        ))
      }))
    }
    diagnostics <- diagnose_chains(chains, warmup = warmup)
  })[["elapsed"]]
  list(
    chains = chains,
    diagnostics = diagnostics,
    warmup = warmup,
    elapsed_sec = elapsed,
    initial_value_seed = initial_value_seed,
    starting_values = starting_values
  )
}

combine_chains <- function(chains, warmup) {
  keep <- lapply(chains, post_warmup_index, warmup = warmup)
  out <- chains[[1]]
  matrix_fields <- c("h.hat", "beta", "lambda", "r", "acc.r", "acc.lambda",
                     "delta", "ystar", "hnew")
  vector_fields <- c("sigsq.eps", "acc.rdelta", "move.type")
  for (field in matrix_fields) {
    if (!is.null(out[[field]])) {
      out[[field]] <- do.call(rbind, lapply(seq_along(chains), function(i) {
        chains[[i]][[field]][keep[[i]], , drop = FALSE]
      }))
    }
  }
  for (field in vector_fields) {
    if (!is.null(out[[field]])) {
      out[[field]] <- unlist(lapply(seq_along(chains), function(i) chains[[i]][[field]][keep[[i]]]))
    }
  }
  out$iter <- sum(vapply(keep, length, integer(1)))
  out
}

safe_real <- function(expr) {
  tryCatch(as.numeric(expr), error = function(...) NA_real_)
}

overall_q25_to_q75 <- function(fit) {
  safe_real({
    overall <- bkmr::OverallRiskSummaries(
      fit = fit,
      qs = c(0.25, 0.75),
      q.fixed = 0.5,
      method = "exact",
      sel = seq_len(fit$iter)
    )
    high <- overall$est[overall$quantile == 0.75]
    low <- overall$est[overall$quantile == 0.25]
    high[1] - low[1]
  })
}

single_var_estimate <- function(fit, exposure = "z1") {
  idx <- match(exposure, colnames(fit$Z))
  if (is.na(idx)) {
    return(NA_real_)
  }
  safe_real({
    uni <- bkmr::SingVarRiskSummaries(
      fit = fit,
      which.z = idx,
      q.fixed = 0.5,
      method = "exact",
      sel = seq_len(fit$iter)
    )
    as.numeric(uni$est[1])
  })
}

top_pip <- function(fit) {
  if (!isTRUE(fit$varsel)) {
    return(data.frame(top_pip_var = NA_character_, top_pip = NA_real_))
  }
  out <- tryCatch(bkmr::ExtractPIPs(fit), error = function(...) NULL)
  if (is.null(out) || nrow(out) == 0L) {
    return(data.frame(top_pip_var = NA_character_, top_pip = NA_real_))
  }
  out <- out[order(-out$PIP), , drop = FALSE]
  data.frame(top_pip_var = as.character(out$variable[1]), top_pip = as.numeric(out$PIP[1]))
}

diagnostic_summary <- function(diag_obj) {
  tab <- diag_obj$table
  rhat <- tab$rhat[is.finite(tab$rhat)]
  bulk_ess <- tab$bulk_ess[is.finite(tab$bulk_ess)]
  tail_ess <- tab$tail_ess[is.finite(tab$tail_ess)]
  min_ess <- pmin(tab$bulk_ess, tab$tail_ess, na.rm = TRUE)
  min_ess <- min_ess[is.finite(min_ess)]
  mcse <- tab$mcse_mean[is.finite(tab$mcse_mean)]
  data.frame(
    max_rhat = if (length(rhat)) max(rhat) else NA_real_,
    max_rank_normalized_rhat = if (length(rhat)) max(rhat) else NA_real_,
    min_ess = if (length(min_ess)) min(min_ess) else NA_real_,
    min_bulk_ess = if (length(bulk_ess)) min(bulk_ess) else NA_real_,
    min_tail_ess = if (length(tail_ess)) min(tail_ess) else NA_real_,
    max_mcse = if (length(mcse)) max(mcse) else NA_real_,
    n_rhat_flag = sum(tab$rhat_flag, na.rm = TRUE),
    n_bulk_ess_flag = sum(tab$bulk_ess_flag, na.rm = TRUE),
    n_tail_ess_flag = sum(tab$tail_ess_flag, na.rm = TRUE),
    n_ess_flag = sum(tab$ess_flag, na.rm = TRUE),
    n_acceptance_flag = sum(tab$acceptance_flag, na.rm = TRUE),
    n_hard_flag = sum(tab$hard_flag, na.rm = TRUE),
    n_flagged = sum(tab$needs_attention, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarize_fit <- function(fit_obj, scenario, model = "exact_public", nd = NA_real_) {
  diag_tab <- fit_obj$diagnostics$table
  ess_vals <- pmin(diag_tab$bulk_ess, diag_tab$tail_ess, na.rm = TRUE)
  ess_vals <- ess_vals[is.finite(ess_vals)]
  converged_hard <- !any(diag_tab$hard_flag, na.rm = TRUE)
  tuning_issue <- any(diag_tab$tuning_flag, na.rm = TRUE)
  combined <- combine_chains(fit_obj$chains, fit_obj$warmup)
  cbind(
    data.frame(
      scenario = scenario,
      model = model,
      nd = nd,
      elapsed_sec = fit_obj$elapsed_sec,
      median_ess_per_sec = if (length(ess_vals)) stats::median(ess_vals) / fit_obj$elapsed_sec else NA_real_,
      converged_hard = converged_hard,
      tuning_issue = tuning_issue,
      converged_strict = converged_hard && !tuning_issue,
      converged = converged_hard && !tuning_issue,
      stringsAsFactors = FALSE
    ),
    diagnostic_summary(fit_obj$diagnostics),
    data.frame(
      overall_q25_to_q75 = overall_q25_to_q75(combined),
      prespecified_univar_var = "z1",
      prespecified_univar_est = single_var_estimate(combined, "z1"),
      stringsAsFactors = FALSE
    ),
    top_pip(combined)
  )
}

error_fit_row <- function(scenario, elapsed, error_message, model = "exact_public", nd = NA_real_) {
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
    max_rank_normalized_rhat = NA_real_,
    min_ess = NA_real_,
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    max_mcse = NA_real_,
    n_rhat_flag = NA_integer_,
    n_bulk_ess_flag = NA_integer_,
    n_tail_ess_flag = NA_integer_,
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

fit_task <- function(scenario, dataset_seed, n, hfun, beta_true, iter, warmup,
                     varsel, chain_seeds, keep_parameter_table = TRUE,
                     parameter_dir = NULL, parameter_id = NULL,
                     initial_value_seed = fixed_initial_value_seed()) {
  dat <- make_bkmr_dataset(dataset_seed, n = n, hfun = hfun, beta_true = beta_true)
  initial_value_task_seed <- derive_initial_value_seed(
    initial_value_seed,
    scenario,
    dataset_seed,
    n,
    hfun,
    beta_true,
    iter,
    warmup,
    varsel,
    paste(chain_seeds, collapse = ";")
  )
  fit_obj <- NULL
  elapsed <- system.time({
    fit_obj <- tryCatch(
      fit_chains(
        y = dat$y,
        Z = dat$Z,
        X = dat$X,
        iter = iter,
        warmup = warmup,
        varsel = varsel,
        chain_seeds = chain_seeds,
        initial_value_seed = initial_value_task_seed
      ),
      error = identity
    )
  })[["elapsed"]]

  if (inherits(fit_obj, "error")) {
    row <- error_fit_row(scenario, elapsed, conditionMessage(fit_obj))
  } else {
    row <- summarize_fit(fit_obj, scenario = scenario)
    row$status <- "ok"
    row$error_message <- NA_character_
    if (keep_parameter_table && !is.null(parameter_dir)) {
      dir.create(parameter_dir, recursive = TRUE, showWarnings = FALSE)
      parameter_file <- sprintf(
        "%s__%s__n%d__varsel_%s__iter_%d.csv",
        scenario,
        dataset_seed,
        n,
        tolower(as.character(varsel)),
        iter
      )
      if (!is.null(parameter_id)) {
        parameter_file <- paste(safe_file_component(parameter_id), parameter_file, sep = "__")
      }
      write_csv(fit_obj$diagnostics$table, file.path(parameter_dir, parameter_file))
    }
  }

  row$dataset_seed <- dataset_seed
  row$n <- n
  row$hfun <- hfun
  row$beta_true <- beta_true
  row$iter <- iter
  row$warmup <- warmup
  row$varsel <- varsel
  row$label_sum <- dat$label_sum
  row$chain_seeds <- paste(chain_seeds, collapse = ";")
  row$initial_value_seed <- initial_value_seed
  row$initial_value_task_seed <- initial_value_task_seed
  row
}

summarize_sweep <- function(details, group_cols) {
  groups <- split(details, interaction(details[group_cols], drop = TRUE, lex.order = TRUE))
  out <- lapply(groups, function(df) {
    ok <- df$status == "ok"
    completed <- df[ok, , drop = FALSE]
    rhat <- completed$max_rhat[is.finite(completed$max_rhat)]
    ess <- completed$min_ess[is.finite(completed$min_ess)]
    bulk_ess <- completed$min_bulk_ess[is.finite(completed$min_bulk_ess)]
    tail_ess <- completed$min_tail_ess[is.finite(completed$min_tail_ess)]
    base <- df[1, group_cols, drop = FALSE]
    cbind(
      base,
      data.frame(
        K = nrow(df),
        n_completed = sum(ok),
        crash_rate = mean(!ok),
        hard_failure_rate = if (sum(ok) == 0L) NA_real_ else mean(!completed$converged_hard),
        tuning_issue_rate = if (sum(ok) == 0L) NA_real_ else mean(completed$tuning_issue),
        mean_label_sum = mean(df$label_sum),
        min_label_sum = min(df$label_sum),
        max_label_sum = max(df$label_sum),
        median_max_rhat = if (length(rhat)) stats::median(rhat) else NA_real_,
        worst_max_rhat = if (length(rhat)) max(rhat) else NA_real_,
        median_min_ess = if (length(ess)) stats::median(ess) else NA_real_,
        worst_min_ess = if (length(ess)) min(ess) else NA_real_,
        median_min_bulk_ess = if (length(bulk_ess)) stats::median(bulk_ess) else NA_real_,
        worst_min_bulk_ess = if (length(bulk_ess)) min(bulk_ess) else NA_real_,
        median_min_tail_ess = if (length(tail_ess)) stats::median(tail_ess) else NA_real_,
        worst_min_tail_ess = if (length(tail_ess)) min(tail_ess) else NA_real_,
        n_hard_pass = if (sum(ok) == 0L) 0L else sum(completed$converged_hard),
        stringsAsFactors = FALSE
      )
    )
  })
  do.call(rbind, out)
}

fixed_chain_seeds <- function() {
  c(142709412L, 1107285817L, 1907135512L, 672667230L)
}

seed_matrix <- function(seed, n_sets, chains = 4L) {
  set.seed(seed)
  replicate(n_sets, sample.int(.Machine$integer.max, chains))
}

dataset_seeds <- function(seed, K) {
  set.seed(seed)
  sample.int(.Machine$integer.max, K)
}
