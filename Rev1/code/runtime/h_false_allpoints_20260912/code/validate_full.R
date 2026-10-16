# Public rerun adapter (2026-10-16); historical source in code/original.
options(stringsAsFactors=FALSE)
out<-'h_false_allpoints_20260912'
reg<-read.csv(file.path(out,'registry.csv'))
reg$reused<-tolower(as.character(reg$reused))=='true'
stopifnot(nrow(reg)==200,!anyDuplicated(reg$task_index),all(table(reg$n)==50))
records<-list()
for(i in seq_len(nrow(reg))){
 id<-reg$task_index[i];x<-readRDS(reg$checkpoint[i])
 runner<-reg$runner[i]
 stopifnot(x$task$task_index==id,x$task$n==reg$n[i],!x$est_h,
  identical(tools::md5sum(runner),x$fingerprint$runner),
  identical(tools::md5sum('Main/Code_Data/code/_pure_bkmr_core.R'),x$fingerprint$core),
  identical(tools::md5sum(names(x$fingerprint$packages)),x$fingerprint$packages),
  x$task$iter==10000,x$task$warmup==5000,!x$task$varsel,
  identical(x$h_indices,seq_len(as.integer(x$task$n))))
 oldpath<-sprintf('revision_20260911/checkpoints/task_%03d.rds',id)
 if(x$status=='completed'){
  old<-readRDS(oldpath)
  stopifnot(x$identity_matches,x$exact_core_draw_matches,
   identical(x$starting_values,old$starting_values),identical(x$task,old$task),
   length(x$draws)==11,nrow(x$diagnostics)==6+x$task$n,x$all_h_alignment_verified,
   all(is.finite(as.matrix(x$diagnostics[c('rhat','bulk_ess','tail_ess','mcse_sd_ratio')]))),nrow(x$chain_status)==4,nrow(x$postprocessing)==4,
   all(x$chain_status$status=='completed'))
  for(k in names(old$draws))stopifnot(identical(as.numeric(x$draws[[k]]),as.numeric(old$draws[[k]])))
  for(k in names(x$draws)){
   m<-x$draws[[k]];stopifnot(identical(dim(m),c(5000L,4L)),all(is.finite(m)))
   z<-x$diagnostics[x$diagnostics$parameter==k,]
   q<-quantile(as.vector(m),c(.025,.975),names=FALSE)
   stopifnot(abs(mean(m)-z$mean)<1e-12,abs(sd(as.vector(m))-z$sd)<1e-12,
    max(abs(q-c(z$q025,z$q975)))<1e-12)
   if(grepl('^h\\[',k))stopifnot(identical(rownames(m),paste0('iter',5001:10000)))
  }
  core<-x$diagnostics[!grepl('^h\\[',x$diagnostics$parameter),]
  primary<-all(is.finite(as.matrix(core[c('rhat','bulk_ess','tail_ess')])))&&all(core$rhat<=1.01&core$bulk_ess>=400&core$tail_ess>=400)
  stopifnot(identical(primary,x$primary_pass))
  a<-read.csv(file.path(reg$source[i],'audit',sprintf('task_%03d_esth_FALSE_identity.csv',id)))
  stopifnot(nrow(a)==30,all(a$matches))
  records[[i]]<-data.frame(task_index=id,n=x$task$n,status=x$status,reused=reg$reused[i],
   primary_pass=primary,exact_core_draws=TRUE,exact_starts=TRUE,h_series=as.integer(x$task$n),
   iterations_per_chain=5000L,chains=4L,aligned_h_iterations=TRUE,
   first_draw_check_max_error=max(x$postprocessing$first_draw_max_abs_error),
   verification_method=if(reg$reused[i])'independent direct draw' else 'independent moments and exact package draw replay',
   independent_mean_max_error=if(reg$reused[i])NA_real_ else max(x$postprocessing$independent_mean_max_abs_error),
   independent_covariance_max_error=if(reg$reused[i])NA_real_ else max(x$postprocessing$independent_covariance_max_abs_error),
   sampling_seconds=sum(x$chain_status$elapsed_seconds),postprocess_seconds=sum(x$postprocessing$postprocess_seconds),
   error_chain=NA_integer_,error_iteration=NA_integer_,error_message='')
 }else{
  stopifnot(id==174,x$original_status=='error',nrow(x$chain_status)==4,
   x$error_context$chain==4,x$error_context$iteration==165,
   x$error_context$message=='H must be positive definite')
  records[[i]]<-data.frame(task_index=id,n=x$task$n,status=x$status,reused=reg$reused[i],
   primary_pass=NA,exact_core_draws=NA,exact_starts=NA,h_series=0L,
   iterations_per_chain=NA_integer_,chains=NA_integer_,aligned_h_iterations=NA,
   first_draw_check_max_error=NA_real_,verification_method='not evaluable',independent_mean_max_error=NA_real_,independent_covariance_max_error=NA_real_,sampling_seconds=sum(x$chain_status$elapsed_seconds),
   postprocess_seconds=NA_real_,error_chain=4L,error_iteration=165L,error_message=x$error_context$message)
 }
}
v<-do.call(rbind,records)
stopifnot(sum(v$status=='completed')==199,sum(v$h_series)==29680,sum(v$primary_pass,na.rm=TRUE)==22)
write.csv(v,file.path(out,'audit/validation.csv'),row.names=FALSE)
cat('PASS: rerun checkpoints validated against freshly generated baseline draws and published core diagnostics.\n')
