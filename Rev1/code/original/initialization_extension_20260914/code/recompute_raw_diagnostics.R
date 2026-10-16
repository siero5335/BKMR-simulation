# Rscript code/recompute_raw_diagnostics.R BUNDLE_ROOT NEW_OUTPUT_DIRECTORY
# A diagnostic audit from raw trajectories; does not rerun the sampler.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop('Usage: Rscript recompute_raw_diagnostics.R BUNDLE_ROOT NEW_OUTPUT_DIRECTORY')
root <- normalizePath(args[1], mustWork=TRUE)
dest <- args[2]
if (dir.exists(dest) || file.exists(dest)) stop('Use a new output directory')
if (as.character(packageVersion('posterior')) != '1.7.0') stop('Use posterior 1.7.0 for this exact audit')
rawdir <- file.path(root, 'raw_chains')
if (!dir.exists(rawdir)) rawdir <- file.path(root, 'chains')
files <- list.files(rawdir, pattern='\\.rds$', full.names=TRUE)
if (!length(files)) stop('No raw trajectories found')
reference <- read.csv(file.path(root, 'parameter_diagnostics.csv'), stringsAsFactors=FALSE)
metrics <- c('mean','sd','q025','q975','rhat','bulk_ess','tail_ess','mcse_mean','mcse_sd_ratio')
meta <- lapply(files, function(path) {
  x <- readRDS(path)
  data.frame(path=path, task_index=x$job$task_index, initial_set=x$job$initial_set,
             chain=x$job$chain, completed=nrow(x$draws))
})
meta <- do.call(rbind,meta)
groups <- unique(meta[c('task_index','initial_set')])
result <- list(); count <- 0L
for (i in seq_len(nrow(groups))) {
  g <- groups[i,]
  mm <- meta[meta$task_index==g$task_index & meta$initial_set==g$initial_set,]
  mm <- mm[order(mm$chain),]
  if (!identical(as.integer(mm$chain),1:4)) stop('A raw group does not contain exactly chains 1:4')
  raw <- lapply(mm$path,readRDS)
  targets <- reference[reference$task_index==g$task_index & reference$initial_set==g$initial_set,]
  for (j in seq_len(nrow(targets))) {
    r <- targets[j,]
    stopifnot(all(mm$completed>=r$budget))
    mat <- do.call(cbind,lapply(raw,function(x)x$draws[seq.int(r$warmup+1L,r$budget),r$parameter]))
    q <- unname(quantile(c(mat),c(.025,.975)))
    se <- posterior::mcse_mean(mat); sdv <- sd(c(mat))
    got <- c(mean=mean(mat),sd=sdv,q025=q[1],q975=q[2],rhat=posterior::rhat(mat),
             bulk_ess=posterior::ess_bulk(mat),tail_ess=posterior::ess_tail(mat),
             mcse_mean=se,mcse_sd_ratio=se/sdv)
    want <- unlist(r[metrics],use.names=TRUE)
    match <- (is.na(got)&is.na(want)) | (is.finite(got)&is.finite(want)&abs(got-want)<=1e-9*pmax(1,abs(want)))
    count <- count+1L
    result[[count]] <- data.frame(task_index=g$task_index,initial_set=g$initial_set,budget=r$budget,
      window=r$window,parameter=r$parameter,metric=metrics,recomputed=as.numeric(got),reference=as.numeric(want),
      absolute_difference=abs(got-want),matches=match)
  }
}
audit <- do.call(rbind,result)
stopifnot(nrow(audit)>0,all(audit$matches))
dir.create(dest,recursive=TRUE)
write.csv(audit,file.path(dest,'raw_metric_comparisons.csv'),row.names=FALSE)
write.csv(data.frame(raw_chains=nrow(meta),initialization_groups=nrow(groups),
   parameter_window_budget_rows=count,metric_comparisons=nrow(audit),all_match=all(audit$matches),
   maximum_absolute_difference=max(audit$absolute_difference,na.rm=TRUE)),
   file.path(dest,'raw_validation_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
cat('Verified',nrow(meta),'raw chains;',nrow(audit),'metric comparisons.\n')
