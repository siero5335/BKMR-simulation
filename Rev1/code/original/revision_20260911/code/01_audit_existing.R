#!/usr/bin/env Rscript
# Run from the BKMR workspace root. All writes are confined to revision_20260911.
options(stringsAsFactors = FALSE)
out <- 'revision_20260911'
src <- 'Main/Code_Data'
source(file.path(src, 'code/_pure_bkmr_core.R'))
wc <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
tasks <- read.csv(file.path(src, 'results/task_table.csv'))
base <- read.csv(file.path(src, 'results/all_task_results.csv'))
stopifnot(!anyDuplicated(tasks$task_index), !anyDuplicated(base$task_index),
          identical(sort(tasks$task_index), sort(base$task_index)))
files <- list.files(file.path(src,'results/parameter_diagnostics_by_task'),
                    pattern='^task_[0-9]+__.*[.]csv$', recursive=TRUE, full.names=TRUE)
target <- c('beta[1]', 'lambda', paste0('r[z',1:4,']'))
pars <- do.call(rbind,lapply(files,function(f) {
  d <- read.csv(f); d <- d[d$parameter %in% target,c('parameter','rhat','bulk_ess','tail_ess','mcse_mean','acceptance')]
  stopifnot(nrow(d)==6L, setequal(d$parameter,target))
  d$task_index <- as.integer(sub('^task_([0-9]+)__.*','\\1',basename(f)))
  d$source_file <- f; d
}))
stopifnot(!anyDuplicated(paste(pars$task_index,pars$parameter)),
          setequal(unique(pars$task_index),base$task_index[base$status=='ok']))
pars$rhat_flag <- !is.finite(pars$rhat) | pars$rhat > 1.01
pars$bulk_flag <- !is.finite(pars$bulk_ess) | pars$bulk_ess < 400
pars$tail_flag <- !is.finite(pars$tail_ess) | pars$tail_ess < 400
pars$flag <- pars$rhat_flag | pars$bulk_flag | pars$tail_flag
wc(pars,'parameter_diagnostics_clean.csv')
flags <- do.call(rbind,lapply(split(pars,pars$task_index),function(d) data.frame(
  task_index=d$task_index[1], beta_flag=d$flag[d$parameter=='beta[1]'],
  lambda_flag=d$flag[d$parameter=='lambda'],any_r_flag=any(d$flag[grepl('^r',d$parameter)]),
  primary_pass=!any(d$flag),rhat_pass=!any(d$rhat_flag),bulk_pass=!any(d$bulk_flag),
  tail_pass=!any(d$tail_flag),both_ess_pass=!any(d$bulk_flag|d$tail_flag))))
flags <- merge(tasks,flags,by='task_index')
wc(flags,'primary_status_by_task.csv')
checks <- data.frame(check=c('prespecified','completed','numerical_noncompletion','primary_pass',
  'rhat_only_pass','bulk_only_pass','tail_only_pass','both_ess_pass','nonfinite_diagnostics','primary_label_disagreements'),
  observed=c(nrow(tasks),nrow(flags),sum(base$status=='error'),sum(flags$primary_pass),
    sum(flags$rhat_pass),sum(flags$bulk_pass),sum(flags$tail_pass),sum(flags$both_ess_pass),
    sum(!is.finite(as.matrix(pars[c('rhat','bulk_ess','tail_ess')]))),
    sum(flags$primary_pass!=base$converged_hard[match(flags$task_index,base$task_index)])),
  expected=c(431,430,1,30,55,85,44,44,0,0))
checks$matches <- checks$observed==checks$expected
wc(checks,'baseline_diagnostics_check.csv')
stopifnot(all(checks$matches))
overlap <- function(d) {
  grid <- expand.grid(any_r_flag=c(FALSE,TRUE),lambda_flag=c(FALSE,TRUE),beta_flag=c(FALSE,TRUE))
  grid <- grid[c('beta_flag','lambda_flag','any_r_flag')]
  grid$tasks <- vapply(seq_len(nrow(grid)),function(i)sum(d$beta_flag==grid$beta_flag[i]&
    d$lambda_flag==grid$lambda_flag[i]&d$any_r_flag==grid$any_r_flag[i]),integer(1)); grid
}
ov <- overlap(flags);wc(ov,'parameter_flag_overlap.csv')
groupcols <- c('bundle','scenario','n','varsel','iter','warmup')
grp <- function(d) split(d, interaction(d[groupcols],drop=TRUE,lex.order=TRUE))
wc(do.call(rbind,lapply(grp(flags),function(d)cbind(d[rep(1,8),groupcols],overlap(d)))),
   'parameter_flag_overlap_by_scenario.csv')
wc(as.data.frame(with(flags,table(beta_flag,lambda_flag))), 'beta_lambda_2x2.csv')
pc <- merge(pars,tasks,by='task_index')
wc(aggregate(cbind(rhat_flag,bulk_flag,tail_flag,flag)~parameter,pc,sum),'parameter_flag_by_metric.csv')
wc(aggregate(cbind(rhat_flag,bulk_flag,tail_flag,flag)~bundle+scenario+n+varsel+iter+parameter,pc,sum),
   'parameter_flag_by_metric_and_scenario.csv')
lam <- pc[pc$parameter=='lambda',c('task_index','dataset_seed',groupcols,'acceptance','rhat')]
wc(lam,'lambda_acceptance_by_task.csv')
corrows <- do.call(rbind,lapply(grp(lam),function(d){
  valid <- is.finite(d$acceptance)&is.finite(d$rhat)
  eligible <- d$iter[1]==10000 && length(unique(d$dataset_seed))>=10 && all(valid) &&
    length(unique(d$acceptance))>1 && length(unique(d$rhat))>1
  cbind(d[1,groupcols],tasks=nrow(d),nonfinite=sum(!valid),
    spearman=if(eligible)cor(d$acceptance,d$rhat,method='spearman') else NA_real_,
    note=if(eligible)'descriptive association across datasets; no p-value' else 'not a repeated-dataset 10000-iteration stratum; correlation not estimated')
}))
wc(corrows,'lambda_acceptance_diagnostics.csv')
png(file.path(out,'lambda_acceptance_scatter.png'),width=1800,height=1500,res=160)
groups <- Filter(function(d)d$iter[1]==10000 && length(unique(d$dataset_seed))>=10,grp(lam))
par(mfrow=c(3,3),mar=c(4,4,3.5,1),mgp=c(2.5,.8,0))
for(d in groups)plot(d$acceptance,d$rhat,pch=19,col='#285A7A99',cex=.65,
  xlab='Mean post-warmup acceptance',ylab='Lambda R-hat',
  main=paste(d$bundle[1],paste0('n=',d$n[1],'  varsel=',d$varsel[1]),
             paste0('iter=',d$iter[1],'  tasks=',nrow(d)),sep='\n'),cex.main=.7)
dev.off()
# Inspect every local RDS, including legacy duplicates, without using their metadata as draws.
rdsfiles <- list.files('old',pattern='[.](rds|RDS|rda|RData|Rdata)$',recursive=TRUE,full.names=TRUE)
inventory <- do.call(rbind,lapply(rdsfiles,function(f){
  x <- tryCatch(readRDS(f),error=identity)
  data.frame(path=f,bytes=file.info(f)$size,object_class=paste(class(x),collapse=';'),
    fields=paste(names(x),collapse=';'),
    has_fit_or_draw_field=any(c('chains','beta','draws','fit','fit_obj') %in% names(x)),
    disposition=if(grepl('_pure_bkmr_parallel_cache',f))'legacy task summary; no draws'
      else 'April numerical tracing; not target June task')
}))
wc(inventory,'audit/rds_inventory.csv')
# Version-specific source snapshots and current environment, no installation.
for(p in c('bkmr','posterior','tmvtnorm'))requireNamespace(p,quietly=TRUE)
capture.output({print(Sys.time());print(sessionInfo());print(RNGkind());
  for(p in c('bkmr','posterior','tmvtnorm'))print(packageDescription(p)[c('Package','Version','RemoteSha','Built')])
  print(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')))
},file=file.path(out,'audit/current_environment.txt'))
for(nm in c('SimData','HFun2','ystar.update.noh','kmbayes','makeVcomps')) {
  f <- tryCatch(getFromNamespace(nm,'bkmr'),error=function(e)NULL)
  if(!is.null(f))capture.output(print(f),file=file.path(out,'audit',paste0('bkmr_',nm,'.R.txt')))
}
for(nm in c('mcse_mean.default','mcse_quantile.default','ess_bulk.default','ess_tail.default'))
  capture.output(print(getFromNamespace(nm,'posterior')),file=file.path(out,'audit',paste0('posterior_',nm,'.R.txt')))
# D1: input generation only, preserving the original data-generating function.
t <- tasks[tasks$task_index==174,]; dat <- make_bkmr_dataset(t$dataset_seed,t$n,t$hfun,t$beta_true)
seed <- derive_initial_value_seed(t$initial_value_seed,t$scenario,t$dataset_seed,t$n,t$hfun,t$beta_true,
                                 t$iter,t$warmup,t$varsel,t$chain_seeds)
stopifnot(dat$label_sum==base$label_sum[base$task_index==174],seed==1740350714L)
saveRDS(list(task=t,initial_value_task_seed=seed,data=dat),file.path(out,'failure_case/input.rds'))
wc(data.frame(row=seq_len(t$n),dat$Z),'failure_case/Z.csv')
wc(data.frame(row=seq_len(t$n),X=dat$X[,1],y=dat$y),'failure_case/X_y.csv')
wc(data.frame(variable=colnames(dat$Z),cor(dat$Z)),'failure_case/exposure_correlation.csv')
inputs <- cbind(dat$Z,X=dat$X[,1],y=dat$y)
wc(data.frame(variable=colnames(inputs),minimum=apply(inputs,2,min),maximum=apply(inputs,2,max),
  sd=apply(inputs,2,sd),nonfinite=colSums(!is.finite(inputs))),'failure_case/input_checks.csv')
input_stats <- function(d){
  s <- svd(scale(d$Z,center=TRUE,scale=FALSE))$d
  data.frame(label_sum=d$label_sum,event_rate=mean(d$y),max_abs_correlation=max(abs(cor(d$Z)[upper.tri(cor(d$Z))])),
    duplicate_rows=sum(duplicated(d$Z)),minimum_pairwise_distance=min(dist(d$Z)),
    centered_rank=qr(scale(d$Z,center=TRUE,scale=FALSE))$rank,
    condition_number=max(s)/min(s),s1=s[1],s2=s[2],s3=s[3],s4=s[4])
}
wc(cbind(task_index=174,initial_value_task_seed=seed,input_stats(dat)),'failure_case/input_geometry.csv')
peers <- base[base$bundle=='balanced_main'&base$n==320&!base$varsel,]
peerstats <- do.call(rbind,lapply(seq_len(nrow(peers)),function(i){
  p <- peers[i,];d<-make_bkmr_dataset(p$dataset_seed,p$n,p$hfun,p$beta_true)
  stopifnot(d$label_sum==p$label_sum)
  cbind(p[c('task_index','dataset_seed','status')],input_stats(d))
}))
wc(peerstats,'failure_case/n320_input_comparison.csv')
writeLines(c('# Input inventory','',
  'Selected baseline: Main/Code_Data. The nominal chain_specific_initial_values_20260621 folder contains prompts only.',
  sprintf('Read %d unique completed-task parameter CSVs; duplicates in other folders were not pooled.',length(files)),
  sprintf('Read all %d local RDS files under old/: no usable target-fit draws found.',nrow(inventory)),
  '862 RDS files contain legacy summary metadata; two are April tracing snapshots for another task.',
  'ZIP member inventories are recorded separately. No target draw files found in archives.',
  'Current and recorded June environment agree in R 4.5.2, macOS 26.3, bkmr 0.2.2.9000 commit 45413e338a316362d629f53bd2a917c4bf485c1e, posterior 1.7.0, tmvtnorm 1.7, BLAS/LAPACK paths.',
  'The June environment file does not record RNGkind or binary hashes; exact pilot diagnostic checks remain necessary.',
  'Secondary dependency versions differ: see audit/environment_comparison.csv. Matching key packages alone does not imply complete environment identity.',
  'A/B require new retained draws from original-condition pilots and a successful identity gate before further tasks.',
  'Original outcomes and source files are retained unchanged.'),file.path(out,'input_inventory.md'))
print(checks);print(ov);print(input_stats(dat))
