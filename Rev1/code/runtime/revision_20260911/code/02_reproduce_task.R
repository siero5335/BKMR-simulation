# Public rerun adapter (2026-10-16); historical source in code/original.
#!/usr/bin/env Rscript
# Usage: Rscript revision_20260911/code/02_reproduce_task.R TASK_INDEX
# Calls the original fitting function without changing its sampler or RNG sequence.
options(stringsAsFactors=FALSE)
out <- 'revision_20260911'; src <- 'Main/Code_Data'
source(file.path(src,'code/_pure_bkmr_core.R'))
id <- as.integer(commandArgs(TRUE)[1]); stopifnot(length(id)==1,!is.na(id))
t <- read.csv(file.path(src,'results/task_table.csv'))
t <- t[t$task_index==id,];stopifnot(nrow(t)==1,t$bundle %in% c('balanced_main','balanced_ablation'),!t$varsel,
  t$n %in% c(40,80,160,320),t$iter==10000,t$warmup==5000,t$hfun==2,t$beta_true==.5)
base <- read.csv(file.path(src,'results/all_task_results.csv')); base<-base[base$task_index==id,]
stopifnot(as.character(packageVersion('posterior'))=='1.7.0',
  packageDescription('bkmr')$RemoteSha=='45413e338a316362d629f53bd2a917c4bf485c1e')
hash_object <- function(x){f<-tempfile();on.exit(unlink(f));saveRDS(x,f,version=3);unname(tools::md5sum(f))}
pkgfiles <- unlist(lapply(c('bkmr','posterior','tmvtnorm','mvtnorm'),function(p)
  list.files(find.package(p),recursive=TRUE,full.names=TRUE)))
fingerprint <- list(task=t,core=tools::md5sum(file.path(src,'code/_pure_bkmr_core.R')),
  runner=tools::md5sum('revision_20260911/code/02_reproduce_task.R'),
  baseline=tools::md5sum(file.path(src,'results/all_task_results.csv')),
  packages=tools::md5sum(pkgfiles),R=R.version,RNGkind=RNGkind(),
  BLAS=extSoftVersion()['BLAS'],LAPACK=La_version(),
  threading=Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')))
hash <- hash_object(fingerprint)
checkpoint <- file.path(out,'checkpoints',sprintf('task_%03d.rds',id))
if(file.exists(checkpoint)){
  old<-readRDS(checkpoint);stopifnot(identical(old$config_hash,hash));
  cat('Verified checkpoint; skipping task',id,'\n');quit(status=0)
}
dat <- make_bkmr_dataset(t$dataset_seed,t$n,t$hfun,t$beta_true)
seeds <- as.integer(strsplit(t$chain_seeds,';',fixed=TRUE)[[1]])
seed <- derive_initial_value_seed(t$initial_value_seed,t$scenario,t$dataset_seed,t$n,t$hfun,
  t$beta_true,t$iter,t$warmup,t$varsel,t$chain_seeds)
stopifnot(dat$label_sum==base$label_sum,seed==base$initial_value_task_seed)
error_context <- NULL
cat('START',format(Sys.time()),'task',id,'n',t$n,'initial seed',seed,'\n');flush.console()
started <- proc.time()[['elapsed']]
fit <- tryCatch(withCallingHandlers(fit_chains(dat$y,dat$Z,dat$X,t$iter,t$warmup,t$varsel,seeds,seed),
  error=function(e){
    # Inspect only stack metadata in the first, untraced execution; no RNG use.
    frames <- sys.frames()
    selected <- lapply(seq_along(frames),function(j){
      env<-frames[[j]];keys<-intersect(c('i','s','iter'),ls(env,all.names=TRUE))
      vals<-lapply(keys,function(k){v<-get(k,env,inherits=FALSE);if(is.atomic(v)&&length(v)<=4)v else NULL})
      names(vals)<-keys;list(frame=j,values=vals)
    })
    error_context <<- list(message=conditionMessage(e),conditionCall=conditionCall(e),
      calls=lapply(sys.calls(),deparse),frame_indices=selected)
  }),error=identity)
elapsed <- proc.time()[['elapsed']]-started
result <- list(task=t,config_hash=hash,fingerprint=fingerprint,initial_value_task_seed=seed,
  label_sum=dat$label_sum,elapsed_sec=elapsed,original_status=base$status,
  current_status=if(inherits(fit,'error'))'error' else 'ok',error_context=error_context)
if(!inherits(fit,'error')){
  d <- diagnostic_draw_table(fit$chains,t$warmup)
  d<-d[d$parameter %in% c('beta[1]','lambda',paste0('r[z',1:4,']')),]
  result$draws <- setNames(d$draws,d$parameter)
  stopifnot(length(result$draws)==6,all(vapply(result$draws,function(x)identical(dim(x),c(5000L,4L)),logical(1))))
  result$starting_values <- fit$starting_values
  result$diagnostics <- fit$diagnostics$table[fit$diagnostics$table$parameter %in% d$parameter,
    c('parameter','rhat','bulk_ess','tail_ess','mcse_mean','acceptance')]
  if(base$status=='ok'){
    origfiles<-list.files(file.path(src,'results/parameter_diagnostics_by_task'),
      pattern=sprintf('^task_%03d__',id),recursive=TRUE,full.names=TRUE)
    stopifnot(length(origfiles)==1)
    original<-read.csv(origfiles);original<-original[match(result$diagnostics$parameter,original$parameter),]
    compare<-do.call(rbind,lapply(c('rhat','bulk_ess','tail_ess','mcse_mean','acceptance'),function(m){
      a<-original[[m]];b<-result$diagnostics[[m]]
      # Only CSV representation-level rounding is tolerated; no statistical tolerance.
      delta<-abs(a-b);tol<-1e-10*pmax(1,abs(a))
      match<-(is.na(a)&is.na(b)) | (is.finite(a)&is.finite(b)&delta<=tol)
      data.frame(task_index=id,parameter=result$diagnostics$parameter,metric=m,original=a,current=b,
        absolute_difference=delta,tolerance=tol,matches=match)
    }))
    result$identity_matches <- all(compare$matches)
    write.csv(compare,file.path(out,'audit',sprintf('task_%03d_identity.csv',id)),row.names=FALSE)
  }
} else result$identity_matches <- FALSE
saveRDS(result,paste0(checkpoint,'.tmp'),version=3)
stopifnot(file.rename(paste0(checkpoint,'.tmp'),checkpoint))
cat('END',format(Sys.time()),'task',id,'status',result$current_status,'elapsed',elapsed,
  'identity',result$identity_matches,'\n')
if(!is.null(error_context))print(error_context)

if(result$current_status=="ok")stopifnot(isTRUE(result$identity_matches))
if(result$current_status=="error")stopifnot(id==174L,result$error_context$message=="H must be positive definite")
