# Public rerun adapter (2026-10-16); historical source in code/original.
options(stringsAsFactors=FALSE)
out<-'initialization_extension_20260914'
reg<-read.csv(file.path(out,'registry.csv'));cases<-read.csv(file.path(out,'cases.csv'))
stopifnot(nrow(reg)==144,nrow(cases)==12,!anyDuplicated(reg$job_id))
chain_status<-list();inputs<-list();loaded<-list()
for(i in seq_len(nrow(reg))){
 r<-reg[i,];x<-readRDS(file.path(out,'chains',paste0(r$job_id,'.rds')))
 stopifnot(identical(x$job,r),identical(x$fingerprint$runner,tools::md5sum(file.path(out,'code/run_chain.R'))),
  identical(x$fingerprint$input,tools::md5sum(file.path(out,'inputs',sprintf('task_%03d.rds',r$task_index)))),
  identical(x$fingerprint$plan,tools::md5sum(file.path(out,'plan.json'))),
  identical(x$fingerprint$packages,tools::md5sum(names(x$fingerprint$packages))),
  x$completed_iterations==nrow(x$draws),ncol(x$draws)==6,all(is.finite(x$draws)))
 inp<-inputs[[as.character(r$task_index)]]
 if(is.null(inp)){inp<-readRDS(file.path(out,'inputs',sprintf('task_%03d.rds',r$task_index)));inputs[[as.character(r$task_index)]]<-inp}
 stopifnot(identical(x$starting_values,inp$initial_sets[[r$initial_set]][[r$chain]]))
 if(r$initial_set=='original'){
  stopifnot(isTRUE(x$original_prefix_match),nrow(x$draws)>=10000)
  for(p in colnames(x$draws))stopifnot(identical(as.numeric(x$draws[5001:10000,p]),as.numeric(inp$reference_draws[[p]][,r$chain])))
 }
 chain_status[[i]]<-data.frame(job_id=r$job_id,task_index=r$task_index,n=r$n,initial_set=r$initial_set,chain=r$chain,
  status=x$status,completed_iterations=x$completed_iterations,error_iteration=x$error_iteration,error_message=x$error_message,
  seconds=x$elapsed_seconds,original_prefix_match=x$original_prefix_match)
 # Loading all raw draws simultaneously is unnecessary; read each four-chain set below.
}
statuses<-do.call(rbind,chain_status);write.csv(statuses,file.path(out,'chain_status.csv'),row.names=FALSE)
metrics<-list();fitrows<-list();chainrows<-list();mi<-0L;fi<-0L;ci<-0L
safe<-function(expr)tryCatch(as.numeric(expr),error=function(e)NA_real_)
passes<-function(d)all(is.finite(as.matrix(d[c('rhat','bulk_ess','tail_ess')])))&&all(d$rhat<=1.01&d$bulk_ess>=400&d$tail_ess>=400)
for(id in cases$task_index)for(set in c('original','alternative_1','alternative_2')){
 g<-reg[reg$task_index==id & reg$initial_set==set,];g<-g[order(g$chain),];n<-g$n[1]
 chs<-lapply(g$job_id,function(k)readRDS(file.path(out,'chains',paste0(k,'.rds'))))
 completed<-vapply(chs,function(x)x$completed_iterations,integer(1))
 later_errors<-sum(vapply(chs,function(x)x$status!='completed',logical(1)))
 for(budget in c(10000L,20000L,40000L,100000L)){
  specs<-list(fixed_5000=5000L)
  if(budget>10000L)specs$half_chain<-as.integer(budget/2)
  for(window in names(specs)){
   warm<-specs[[window]];available<-all(completed>=budget)
   fr<-data.frame(task_index=id,n,initial_set=set,budget,window,warmup=warm,retained_per_chain=budget-warm,
    evaluable=available,chains_with_later_or_current_error=later_errors,minimum_completed_iterations=min(completed),
    primary_pass=NA,beta_pass=NA,beta_mean=NA_real_,beta_sd=NA_real_,beta_q025=NA_real_,beta_q975=NA_real_,
    beta_error=NA_real_,beta_abs_error=NA_real_,beta_contains_truth=NA,beta_contains_zero=NA,
    beta_prob_positive=NA_real_,beta_rhat=NA_real_,beta_bulk_ess=NA_real_,beta_tail_ess=NA_real_,
    beta_mcse=NA_real_,beta_mcse_sd_ratio=NA_real_,chain_mean_range=NA_real_,chain_mean_range_over_sd=NA_real_)
   if(available){
    keep<-seq.int(warm+1L,budget);tab<-list()
    for(p in colnames(chs[[1]]$draws)){
     m<-vapply(chs,function(x)x$draws[keep,p],numeric(length(keep)))
     q<-quantile(as.vector(m),c(.025,.975),names=FALSE);sdv<-sd(as.vector(m));mcse<-safe(posterior::mcse_mean(m))
     z<-data.frame(parameter=p,mean=mean(m),sd=sdv,q025=q[1],q975=q[2],rhat=safe(posterior::rhat(m)),
      bulk_ess=safe(posterior::ess_bulk(m)),tail_ess=safe(posterior::ess_tail(m)),mcse_mean=mcse,mcse_sd_ratio=mcse/sdv)
     tab[[p]]<-z;mi<-mi+1L;metrics[[mi]]<-cbind(task_index=id,n,initial_set=set,budget,window,warmup=warm,z)
     if(p=='beta[1]'){
      cm<-colMeans(m)
      fr$beta_prob_positive<-mean(m>0);fr$chain_mean_range<-diff(range(cm));fr$chain_mean_range_over_sd<-diff(range(cm))/sdv
      for(chain in 1:4){
       cq<-quantile(m[,chain],c(.025,.975),names=FALSE);ci<-ci+1L
       chainrows[[ci]]<-data.frame(task_index=id,n,initial_set=set,chain,budget,window,warmup=warm,
        beta_mean=cm[chain],beta_q025=cq[1],beta_q975=cq[2],beta_error=cm[chain]-.5)
      }
     }
    }
    tt<-do.call(rbind,tab);bb<-tab[['beta[1]']]
    fr$primary_pass<-passes(tt);fr$beta_pass<-passes(bb)
    for(pair in list(c('beta_mean','mean'),c('beta_sd','sd'),c('beta_q025','q025'),c('beta_q975','q975'),
     c('beta_rhat','rhat'),c('beta_bulk_ess','bulk_ess'),c('beta_tail_ess','tail_ess'),c('beta_mcse','mcse_mean'),c('beta_mcse_sd_ratio','mcse_sd_ratio')))fr[[pair[1]]]<-bb[[pair[2]]]
    fr$beta_error<-bb$mean-.5;fr$beta_abs_error<-abs(fr$beta_error)
    fr$beta_contains_truth<-bb$q025<=.5 && bb$q975>=.5;fr$beta_contains_zero<-bb$q025<=0 && bb$q975>=0
   }
   fi<-fi+1L;fitrows[[fi]]<-fr
  }
 }
 rm(chs);invisible(gc(FALSE))
}
fits<-do.call(rbind,fitrows);diag<-do.call(rbind,metrics);chains<-do.call(rbind,chainrows)
stopifnot(nrow(fits)==252)
write.csv(fits,file.path(out,'fit_summaries.csv'),row.names=FALSE)
write.csv(diag,file.path(out,'parameter_diagnostics.csv'),row.names=FALSE)
write.csv(chains,file.path(out,'beta_by_chain.csv'),row.names=FALSE)
# Recheck original 10k diagnostic summaries as well as exact prefix draws.
for(id in cases$task_index){
 old<-inputs[[as.character(id)]]$baseline_diagnostics
 current<-diag[diag$task_index==id & diag$initial_set=='original' & diag$budget==10000 & diag$window=='fixed_5000',]
 current<-current[match(old$parameter,current$parameter),]
 for(k in c('mean','sd','q025','q975','rhat','bulk_ess','tail_ess','mcse_mean'))stopifnot(all(abs(old[[k]]-current[[k]])<=1e-9*pmax(1,abs(old[[k]]))))
}
write.csv(data.frame(chain_jobs=144,original_prefix_checks=48,all_fresh_baseline_prefixes_exact=TRUE,
 original_10k_diagnostic_sets=12,all_original_10k_diagnostics_match=TRUE,
 completed_100k_chains=sum(statuses$status=='completed'),sampler_error_chains=sum(statuses$status!='completed')),
 file.path(out,'audit/validation_summary.csv'),row.names=FALSE)
cat('Validated 144 chain outcomes and 48 freshly regenerated baseline prefixes; evaluated 252 fit/window/budget combinations.\n')
