args<-commandArgs(TRUE)
if(length(args)!=1)stop('Usage: Rscript check_saved_draws.R RERUN_WORKSPACE')
setwd(normalizePath(args[1],mustWork=TRUE))
suppressPackageStartupMessages(library(posterior))
x <- readRDS('initialization_extension_20260914/inputs/task_090.rds')
z <- do.call(rbind,lapply(names(x$reference_draws),function(p){
  d <- x$reference_draws[[p]]
  data.frame(parameter=p,rhat=rhat(d),bulk=ess_bulk(d),tail=ess_tail(d),
             bulk_growth=ess_bulk(d)/ess_bulk(d[1:2500,]),tail_growth=ess_tail(d)/ess_tail(d[1:2500,]),
             mcse_mean=mcse_mean(d),mcse_q025=mcse_quantile(d,probs=.025),mcse_q975=mcse_quantile(d,probs=.975))
}))
print(z,digits=7,row.names=FALSE)
f<-readRDS('revision_20260911/failure_case/error_time_trace.rds')$snapshot
str(f$previous_state)
cat('underflow bounds',expm1(1075*log(2)/c(160,320,640)),'\n')
if(!is.null(f$previous_state$lambda))cat('upper bound',320*log1p(f$previous_state$lambda),'\n')
