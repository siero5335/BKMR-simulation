#!/usr/bin/env Rscript
# Refuses to combine original labels with mismatched or incomplete reproduced draws.
options(stringsAsFactors=FALSE)
out <- 'revision_20260911'
stopifnot(as.character(packageVersion('posterior'))=='1.7.0')
tasks <- read.csv('Main/Code_Data/results/task_table.csv')
base <- read.csv('Main/Code_Data/results/all_task_results.csv')
flags <- read.csv(file.path(out,'primary_status_by_task.csv'))
tasks <- subset(tasks,bundle=='balanced_main' & !varsel & n %in% c(80,160,320) & iter==10000 & warmup==5000)
status <- do.call(rbind,lapply(tasks$task_index,function(id){
  p<-file.path(out,'checkpoints',sprintf('task_%03d.rds',id)); x<-if(file.exists(p))readRDS(p) else NULL
  data.frame(task_index=id,n=tasks$n[tasks$task_index==id],original_status=base$status[base$task_index==id],
    checkpoint_exists=!is.null(x),current_status=if(is.null(x))'not_run' else x$current_status,
    identity_matches=if(base$status[base$task_index==id]!='ok'||is.null(x)||is.null(x$identity_matches))NA else x$identity_matches)
}))
write.csv(status,file.path(out,'beta_acquisition_status.csv'),row.names=FALSE)
needed<-status[status$original_status=='ok',]
if(!all(needed$checkpoint_exists & needed$current_status=='ok' & !is.na(needed$identity_matches) & needed$identity_matches))
  stop('A/B aggregation withheld: every original completed task must have identity-verified draws. See beta_acquisition_status.csv.')
rows <- lapply(needed$task_index,function(id){
  x<-readRDS(file.path(out,'checkpoints',sprintf('task_%03d.rds',id)))
  stopifnot(identical(tools::md5sum(names(x$fingerprint$packages)),x$fingerprint$packages))
  m<-x$draws[['beta[1]']]; stopifnot(identical(dim(m),c(5000L,4L)),all(is.finite(m)))
  qs<-quantile(as.vector(m),c(.025,.975),type=7,names=FALSE)
  s<-sd(as.vector(m));mu<-mean(m)
  qmcse<-as.numeric(posterior::mcse_quantile(m,probs=c(.025,.975)))
  mcse<-as.numeric(posterior::mcse_mean(m))
  d<-x$diagnostics[x$diagnostics$parameter=='beta[1]',]
  primary<-flags$primary_pass[match(id,flags$task_index)]
  data.frame(task_index=id,n=x$task$n,primary_status=if(primary)'primary pass' else 'primary non-pass',
    beta_only_pass=all(is.finite(c(d$rhat,d$bulk_ess,d$tail_ess)))&&d$rhat<=1.01&&d$bulk_ess>=400&&d$tail_ess>=400,
    posterior_mean=mu,posterior_sd=s,q025=qs[1],q975=qs[2],interval_width=diff(qs),
    estimation_error=mu-.5,covered=qs[1]<=.5 && .5<=qs[2],
    rhat=posterior::rhat(m),bulk_ess=posterior::ess_bulk(m),tail_ess=posterior::ess_tail(m),
    mcse_mean=mcse,mcse_mean_over_sd=if(s>0)mcse/s else NA_real_,
    mcse_q025=qmcse[1],mcse_q975=qmcse[2],
    mcse_note=if(s<=0)'zero posterior SD' else if(any(!is.finite(c(mcse,qmcse))))
      'nonfinite official MCSE estimate; retained in counts' else 'finite estimate; not a precision guarantee')
})
d<-do.call(rbind,rows)
wc<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE)
wc(d,'beta_estimand_by_task.csv')
wc(d[c('task_index','n','primary_status','beta_only_pass','rhat','bulk_ess','tail_ess',
       'posterior_sd','mcse_mean','mcse_mean_over_sd','mcse_q025','mcse_q975','mcse_note')],'beta_mcse_by_task.csv')
expanded<-rbind(d,transform(d,primary_status='all completed'))
groups<-split(expanded,interaction(expanded$n,expanded$primary_status,drop=TRUE))
summaries<-do.call(rbind,lapply(groups,function(g){
  m<-nrow(g);k<-sum(g$covered);n<-g$n[1];status<-g$primary_status[1]
  z<-qnorm(.975);den<-1+z*z/m;center<-(k/m+z*z/(2*m))/den
  half<-z*sqrt((k/m)*(1-k/m)/m+z*z/(4*m*m))/den
  data.frame(n=n,primary_status=status,prespecified_at_n=50,completed_at_n=sum(d$n==n),
    evaluable_in_stratum=m,numerical_noncompletion_at_n=sum(base$bundle=='balanced_main'&base$n==n&base$status=='error'),
    bias=mean(g$estimation_error),simulation_bias_mcse=if(m>=2)sd(g$estimation_error)/sqrt(m) else NA_real_,
    RMSE=sqrt(mean(g$estimation_error^2)),coverage_k=k,coverage_m=m,coverage_percent=100*k/m,
    simulation_coverage_wilson_lower=center-half,simulation_coverage_wilson_upper=center+half,
    mean_interval_width=mean(g$interval_width))
}))
wc(summaries,'beta_estimand_by_n_and_primary_status.csv')
ms<-do.call(rbind,lapply(groups,function(g)do.call(rbind,lapply(
  c('mcse_mean_over_sd','mcse_q025','mcse_q975'),function(v){
    x<-g[[v]];ok<-is.finite(x)
    data.frame(n=g$n[1],primary_status=g$primary_status[1],metric=v,tasks=nrow(g),
      finite=sum(ok),nonfinite=sum(!ok),median=if(all(ok))median(x) else NA_real_,
      minimum=if(all(ok))min(x) else NA_real_,maximum=if(all(ok))max(x) else NA_real_,
      note=if(all(ok))'all tasks included' else 'summary withheld because some estimates are nonfinite')
  }))))
wc(ms,'beta_mcse_summary.csv')
print(summaries)
