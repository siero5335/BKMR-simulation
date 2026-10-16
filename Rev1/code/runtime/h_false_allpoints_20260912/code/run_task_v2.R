# Public rerun adapter (2026-10-16); historical source in code/original.
options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE);id<-as.integer(args[1]);est_h<-FALSE
out<-'h_false_allpoints_20260912';source('Main/Code_Data/code/_pure_bkmr_core.R')
t<-read.csv('Main/Code_Data/results/task_table.csv');t<-t[t$task_index==id,]
stopifnot(nrow(t)==1,t$bundle %in% c('balanced_main','balanced_ablation'),!t$varsel,
 t$n %in% c(40,80,160,320),t$iter==10000,t$warmup==5000,t$hfun==2,t$beta_true==.5,
 as.character(packageVersion('posterior'))=='1.7.0',
 packageDescription('bkmr')$RemoteSha=='45413e338a316362d629f53bd2a917c4bf485c1e')
reference<-readRDS(sprintf('revision_20260911/checkpoints/task_%03d.rds',id))
stopifnot(identical(tools::md5sum(names(reference$fingerprint$packages)),reference$fingerprint$packages))
stopifnot(id %in% 51:250)
fingerprint<-list(task=t,est_h=est_h,runner=tools::md5sum('h_false_allpoints_20260912/code/run_task_v2.R'),
 core=tools::md5sum('Main/Code_Data/code/_pure_bkmr_core.R'),packages=reference$fingerprint$packages,
 R=R.version,RNGkind=RNGkind(),BLAS=extSoftVersion()['BLAS'],LAPACK=La_version())
prefix<-sprintf('task_%03d_esth_%s',id,est_h)
path<-file.path(out,'checkpoints',paste0(prefix,'.rds'))
if(file.exists(path)){
 old<-readRDS(path);stopifnot(identical(old$fingerprint,fingerprint));cat('SKIP',prefix,'\n');quit(status=0)
}
dat<-make_bkmr_dataset(t$dataset_seed,t$n,t$hfun,t$beta_true)
init<-derive_initial_value_seed(t$initial_value_seed,t$scenario,t$dataset_seed,t$n,t$hfun,t$beta_true,t$iter,t$warmup,t$varsel,t$chain_seeds)
starts<-make_bkmr_binomial_starting_values(dat$y,dat$Z,dat$X,4,FALSE,init)
seeds<-as.integer(strsplit(t$chain_seeds,';',fixed=TRUE)[[1]])
base<-read.csv('Main/Code_Data/results/all_task_results.csv');base<-base[base$task_index==id,]
stopifnot(init==base$initial_value_task_seed,dat$label_sum==base$label_sum)
oldpath<-sprintf('revision_20260911/checkpoints/task_%03d.rds',id)
if(file.exists(oldpath)&&base$status=='ok')stopifnot(identical(starts,readRDS(oldpath)$starting_values))
h_indices<-seq_len(as.integer(t$n))
retained_h_indices<-as.integer(c(1,t$n*c(.25,.5,.75),t$n))
fits<-list();statuses<-list();post_status<-list();error_context<-NULL;start_total<-proc.time()[['elapsed']]
for(i in 1:4){
 set.seed(seeds[i]);started<-proc.time()[['elapsed']];error_iteration<-NA_integer_
 fit<-tryCatch(withCallingHandlers({
  ans<-NULL;invisible(capture.output(ans<-suppressMessages(bkmr::kmbayes(y=dat$y,Z=dat$Z,X=dat$X,
   iter=t$iter,family='binomial',varsel=FALSE,verbose=FALSE,est.h=est_h,starting.values=starts[[i]]))))
  ans
 },error=function(e){
  for(fr in sys.frames())if(exists('s',envir=fr,inherits=FALSE)){
   v<-get('s',envir=fr,inherits=FALSE);if(is.numeric(v)&&length(v)==1)error_iteration<<-as.integer(v)
  }
  error_context<<-list(message=conditionMessage(e),call=conditionCall(e),chain=i,iteration=error_iteration)
 }),error=identity)
 ok<-!inherits(fit,'error')
 statuses[[i]]<-data.frame(task_index=id,n=t$n,est_h,chain=i,status=if(ok)'completed' else 'error',
  elapsed_seconds=proc.time()[['elapsed']]-started,error_iteration,error_message=if(ok)'' else conditionMessage(fit))
 write.csv(do.call(rbind,statuses),file.path(out,'logs',paste0(prefix,'_chains.csv')),row.names=FALSE)
 if(!ok)break
 stopifnot(all(is.finite(fit$beta)),all(is.finite(fit$lambda)),all(is.finite(fit$r)),all(is.finite(fit$ystar)))
 # The sampler has finished; postprocessing never feeds h back into it.
 rng_after_fit<-.Random.seed
 post_seed<-as.integer(900000+10*id+i)
 post_started<-proc.time()[['elapsed']]
 set.seed(post_seed)
 hs<-bkmr::SamplePred(fit,Znew=dat$Z[h_indices,,drop=FALSE],
  Xnew=matrix(0,nrow=t$n,ncol=ncol(dat$X)),sel=5001:10000,type='link')
 stopifnot(identical(dim(hs),c(5000L,as.integer(t$n))),all(is.finite(hs)))
 # Independently reconstruct the first conditional Gaussian draw, using solve()
 # and distances instead of the package's covariance helpers.
 s<-5001L;rr<-as.numeric(fit$r[s,]);lam<-as.numeric(fit$lambda[s,])[1]
 weighted<-sweep(dat$Z,2,sqrt(rr),'*')
 K<-exp(-as.matrix(dist(weighted))^2)
 V<-diag(nrow(K))+lam*K
 resid<-fit$ystar[s,]-drop(dat$X %*% fit$beta[s,])
 Ksub<-K[h_indices,,drop=FALSE]
 mu<-drop(lam*Ksub %*% solve(V,resid))
 Sigma<-as.numeric(fit$sigsq.eps[s])*lam*(K[h_indices,h_indices]-lam*Ksub %*% solve(V,t(Ksub)))
 # Near-singular covariance matrices amplify tiny rounding differences in their roots.
 # Independently verify conditional moments, then exactly replay the package draw.
 K1<-exp(-bkmr:::makeKpart(rr,dat$Z))
 K10<-exp(-bkmr:::makeKpart(rr,dat$Z,dat$Z))
 vv<-bkmr:::makeVcomps(r=rr,lambda=lam,Z=dat$Z,data.comps=fit$data.comps)$Vinv
 prod<-lam*K10 %*% vv
 pm<-prod %*% resid
 pc<-lam*as.numeric(fit$sigsq.eps[s])*(K1-prod %*% t(K10))
 meanerr<-max(abs(mu-pm));coverr<-max(abs(Sigma-pc))
 stopifnot(meanerr<=1e-8*max(1,abs(mu)),coverr<=1e-8*max(1,abs(Sigma)))
 proot<-try(chol(pc),silent=TRUE)
 pmethod<-if(inherits(proot,'try-error'))'svd' else 'chol'
 if(pmethod=='svd'){ss<-svd(pc);proot<-t(ss$v %*% (t(ss$u)*sqrt(ss$d)))}
 set.seed(post_seed);replay<-pm+crossprod(proot,rnorm(t$n))
 maxerr<-max(abs(replay-hs[1,]))
 stopifnot(maxerr<=1e-12*max(1,abs(hs[1,])))
 .Random.seed<-rng_after_fit
 post_status[[i]]<-data.frame(task_index=id,n=t$n,chain=i,seed=post_seed,
  postprocess_seconds=proc.time()[['elapsed']]-post_started,first_draw_max_abs_error=maxerr,independent_mean_max_abs_error=meanerr,independent_covariance_max_abs_error=coverr,first_draw_root_method=pmethod)
 write.csv(do.call(rbind,post_status),file.path(out,'audit',paste0(prefix,'_postprocess.csv')),row.names=FALSE)
 # All h draws remain in memory until their chain-specific diagnostics are computed.
 fits[[i]]<-list(beta=fit$beta,lambda=fit$lambda,r=fit$r,iter=fit$iter,Z=fit$Z,X=fit$X,
  h_selected=hs,acc.lambda=fit$acc.lambda,acc.r=fit$acc.r)
 rm(fit,hs)
 cat(prefix,'chain',i,'completed\n');flush.console()
}
result<-list(task=t,est_h=est_h,fingerprint=fingerprint,initial_value_task_seed=init,
 label_sum=dat$label_sum,starting_values=starts,h_indices=h_indices,retained_h_indices=retained_h_indices,original_status=base$status,
 status=if(length(fits)==4)'completed' else 'error',error_context=error_context,
 chain_status=do.call(rbind,statuses),elapsed_seconds=proc.time()[['elapsed']]-start_total)
if(result$status=='completed'){
 d<-diagnostic_draw_table(fits,t$warmup)
 d<-d[d$parameter %in% c('beta[1]','lambda',paste0('r[z',1:4,']')),]
 for(j in seq_along(h_indices))d<-rbind(d,data.frame(parameter=paste0('h[',h_indices[j],']'),
  draws=I(list(sapply(fits,function(f)f$h_selected[,j])))))
 result$draws<-setNames(d$draws,d$parameter)
 old_draw_path<-sprintf('revision_20260911/checkpoints/task_%03d.rds',id)
 old_draws<-readRDS(old_draw_path)$draws
 result$exact_core_draw_matches<-all(vapply(names(old_draws),function(k)
  isTRUE(all.equal(result$draws[[k]],old_draws[[k]],tolerance=0)),logical(1)))
 stopifnot(result$exact_core_draw_matches)
 result$postprocessing<-do.call(rbind,post_status)
 tab<-do.call(rbind,lapply(seq_len(nrow(d)),function(j){
  m<-d$draws[[j]];q<-quantile(as.vector(m),c(.025,.975),names=FALSE)
  data.frame(parameter=d$parameter[j],mean=mean(m),sd=sd(as.vector(m)),q025=q[1],q975=q[2],
   rhat=posterior::rhat(m),bulk_ess=posterior::ess_bulk(m),tail_ess=posterior::ess_tail(m),
   mcse_mean=posterior::mcse_mean(m),mcse_sd_ratio=posterior::mcse_mean(m)/sd(as.vector(m)))
 }))
 ac<-acceptance_table(fits,t$warmup);tab$acceptance<-ac$acceptance[match(tab$parameter,ac$parameter)]
 stopifnot(all(is.finite(as.matrix(tab[c('mean','sd','q025','q975','rhat','bulk_ess','tail_ess','mcse_mean','mcse_sd_ratio')]))))
 result$all_h_alignment_verified<-all(vapply(result$draws[paste0('h[',h_indices,']')],function(m) identical(dim(m),c(5000L,4L)) && identical(rownames(m),paste0('iter',5001:10000)),logical(1)))
 stopifnot(result$all_h_alignment_verified)
 result$diagnostics<-tab
 core<-tab[!grepl('^h\\[',tab$parameter),]
 result$primary_pass<-all(is.finite(as.matrix(core[c('rhat','bulk_ess','tail_ess')])))&&
  all(core$rhat<=1.01 & core$bulk_ess>=400 & core$tail_ess>=400)
 m<-result$draws[['beta[1]']];result$beta_quantile_mcse<-as.numeric(posterior::mcse_quantile(m,probs=c(.025,.975)))
 if(!est_h){
  original_path<-list.files('Main/Code_Data/results/parameter_diagnostics_by_task',pattern=sprintf('^task_%03d__',id),recursive=TRUE,full.names=TRUE)
  stopifnot(length(original_path)==1)
  original<-read.csv(original_path);original<-original[match(core$parameter,original$parameter),]
  comparisons<-do.call(rbind,lapply(c('rhat','bulk_ess','tail_ess','mcse_mean','acceptance'),function(metric){
   a<-original[[metric]];b<-core[[metric]];tol<-1e-10*pmax(1,abs(a))
   data.frame(parameter=core$parameter,metric,original=a,current=b,matches=(is.na(a)&is.na(b))|(is.finite(a)&is.finite(b)&abs(a-b)<=tol))
  }))
  result$identity_matches<-all(comparisons$matches)
  write.csv(comparisons,file.path(out,'audit',paste0(prefix,'_identity.csv')),row.names=FALSE)
  stopifnot(result$identity_matches)
 }
 write.csv(tab,file.path(out,'audit',paste0(prefix,'_diagnostics.csv')),row.names=FALSE)
}
# Keep a prespecified subset for independent final audit; discard bulky all-point histories.
if(result$status=='completed') result$draws<-result$draws[c('beta[1]','lambda',paste0('r[z',1:4,']'),paste0('h[',retained_h_indices,']'))]
result$elapsed_seconds<-proc.time()[['elapsed']]-start_total
saveRDS(result,paste0(path,'.tmp'));stopifnot(file.rename(paste0(path,'.tmp'),path))
write.csv(data.frame(task_index=id,n=t$n,est_h,status=result$status,
 primary_pass=if(result$status=='completed')result$primary_pass else NA,
 elapsed_seconds=result$elapsed_seconds,error_message=if(is.null(error_context))'' else error_context$message),
 file.path(out,'audit',paste0(prefix,'_status.csv')),row.names=FALSE)
cat(prefix,'DONE',result$status,'elapsed',result$elapsed_seconds,'\n')

if(result$status!='completed')stopifnot(id==174L,result$error_context$chain==4L,result$error_context$iteration==165L,result$error_context$message=='H must be positive definite')
