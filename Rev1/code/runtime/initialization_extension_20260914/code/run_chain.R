options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE);stopifnot(length(args)==1)
out<-'initialization_extension_20260914';job_id<-args[1]
registry<-read.csv(file.path(out,'registry.csv'));job<-registry[registry$job_id==job_id,]
stopifnot(nrow(job)==1)
input_path<-file.path(out,'inputs',sprintf('task_%03d.rds',job$task_index))
inp<-readRDS(input_path);dat<-inp$data;total<-100000L
stopifnot(as.character(packageVersion('posterior'))=='1.7.0',
 packageDescription('bkmr')$RemoteSha=='45413e338a316362d629f53bd2a917c4bf485c1e',
 identical(tools::md5sum(names(inp$fingerprint$packages)),inp$fingerprint$packages))
do.call(RNGkind,as.list(inp$fingerprint$RNGkind))
fingerprint<-list(runner=tools::md5sum(file.path(out,'code/run_chain.R')),input=tools::md5sum(input_path),
 plan=tools::md5sum(file.path(out,'plan.json')),packages=inp$fingerprint$packages,R=R.version,
 RNGkind=RNGkind(),BLAS=extSoftVersion()['BLAS'],iterations=total,job=job)
dest<-file.path(out,'chains',paste0(job_id,'.rds'))
if(file.exists(dest)){
 old<-readRDS(dest);stopifnot(identical(old$fingerprint,fingerprint));cat('SKIP',job_id,'\n');quit(status=0)
}
extract_core<-function(fit,last){
 if(last<1L)return(matrix(numeric(),nrow=0L,ncol=6L,dimnames=list(NULL,c('beta[1]','lambda',paste0('r[z',1:4,']')))))
 keep<-seq_len(last)
 ans<-cbind(as.matrix(fit$beta)[keep,1],as.matrix(fit$lambda)[keep,1],as.matrix(fit$r)[keep,1:4,drop=FALSE])
 colnames(ans)<-c('beta[1]','lambda',paste0('r[z',1:4,']'));ans
}
set.seed(job$sampler_seed)
started<-proc.time()[['elapsed']];error_iteration<-NA_integer_;recovered<-NULL;recovery_problem<-''
fit<-tryCatch(withCallingHandlers({
 ans<-NULL;invisible(capture.output(ans<-suppressMessages(bkmr::kmbayes(y=dat$y,Z=dat$Z,X=dat$X,
  iter=total,family='binomial',varsel=FALSE,est.h=FALSE,verbose=FALSE,
  starting.values=inp$initial_sets[[job$initial_set]][[job$chain]]))))
 ans
},error=function(e){
 for(fr in sys.frames()){
  if(exists('chain',fr,inherits=FALSE)&&exists('s',fr,inherits=FALSE)&&exists('nsamp',fr,inherits=FALSE)){
   s<-get('s',fr,inherits=FALSE);ch<-get('chain',fr,inherits=FALSE)
   if(is.numeric(s)&&length(s)==1&&is.list(ch)&&all(c('beta','lambda','r') %in% names(ch))){
    error_iteration<<-as.integer(s)
    rr<-try(extract_core(ch,max(0L,as.integer(s)-1L)),silent=TRUE)
    if(inherits(rr,'try-error'))recovery_problem<<-as.character(rr) else recovered<<-rr
   }
  }
 }
}),error=identity)
elapsed<-proc.time()[['elapsed']]-started
ok<-!inherits(fit,'error')
if(ok){draws<-extract_core(fit,total);message_text<-''}else{
 draws<-if(is.null(recovered))extract_core(NULL,0L) else recovered
 message_text<-conditionMessage(fit)
}
stopifnot(all(is.finite(draws)),ncol(draws)==6L)
prefix_match<-NA
if(job$initial_set=='original' && nrow(draws)>=10000){
 prefix_match<-all(vapply(colnames(draws),function(k)
  identical(as.numeric(draws[5001:10000,k]),as.numeric(inp$reference_draws[[k]][,job$chain])),logical(1)))
 stopifnot(prefix_match)
}
if(job$initial_set=='original')stopifnot(nrow(draws)>=10000)
result<-list(job=job,fingerprint=fingerprint,status=if(ok)'completed' else 'sampler_error',
 completed_iterations=nrow(draws),error_iteration=error_iteration,error_message=message_text,
 recovery_problem=recovery_problem,elapsed_seconds=elapsed,original_prefix_match=prefix_match,
 starting_values=inp$initial_sets[[job$initial_set]][[job$chain]],draws=draws)
saveRDS(result,paste0(dest,'.tmp'));stopifnot(file.rename(paste0(dest,'.tmp'),dest))
write.csv(data.frame(job_id,task_index=job$task_index,n=job$n,initial_set=job$initial_set,chain=job$chain,
 status=result$status,completed_iterations=nrow(draws),error_iteration,error_message=message_text,
 elapsed_seconds=elapsed,original_prefix_match=prefix_match),file.path(out,'audit',paste0(job_id,'_status.csv')),row.names=FALSE)
cat(job_id,result$status,'iterations',nrow(draws),'seconds',elapsed,'\n')
