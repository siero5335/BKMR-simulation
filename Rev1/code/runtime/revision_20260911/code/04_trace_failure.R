#!/usr/bin/env Rscript
# Separate process, allowed only after the original untraced error is reproduced.
# Error-time stack inspection consumes no random numbers and changes no sampler input.
out<-'revision_20260911';source('Main/Code_Data/code/_pure_bkmr_core.R')
original<-readRDS(file.path(out,'checkpoints/task_174.rds'))
stopifnot(original$current_status=='error',
          original$error_context$message=='H must be positive definite')
input<-readRDS(file.path(out,'failure_case/input.rds'));t<-input$task;dat<-input$data
seeds<-as.integer(strsplit(t$chain_seeds,';',fixed=TRUE)[[1]])
snapshot<-NULL;inspection_error<-NULL
result<-tryCatch(withCallingHandlers(
  fit_chains(dat$y,dat$Z,dat$X,t$iter,t$warmup,t$varsel,seeds,input$initial_value_task_seed),
  error=function(e){
    rng<-get('.Random.seed',envir=.GlobalEnv)
    frames<-sys.frames();calls<-sys.calls()
    tryCatch({
      has<-function(env,keys)all(vapply(keys,exists,logical(1),envir=env,inherits=FALSE))
      hidx<-which(vapply(frames,has,logical(1),keys=c('H','mean','algorithm')))
      kidx<-which(vapply(frames,has,logical(1),keys=c('chain','s','Vcomps','Z')))
      fidx<-which(vapply(frames,has,logical(1),keys=c('chains','chain_seeds','starting_values','i')))
      snapshot<<-list(message=conditionMessage(e),conditionCall=conditionCall(e),
        calls=lapply(calls,deparse),rng=rng,task=t,environment=sessionInfo())
      if(length(hidx)){
        env<-frames[[tail(hidx,1)]]
        snapshot$rtmvnorm_arguments<<-mget(c('H','mean','lower','upper','algorithm'),env,inherits=FALSE)
        snapshot$rtmvnorm_arguments$start.value<<-eval(quote(list(...)),env)$start.value
      }
      if(length(kidx)){
        env<-frames[[tail(kidx,1)]];s<-get('s',env);ch<-get('chain',env)
        snapshot$iteration<<-s
        snapshot$previous_state<<-list(lambda=ch$lambda[s-1,,drop=FALSE],r=ch$r[s-1,,drop=FALSE],
          beta=ch$beta[s-1,,drop=FALSE],ystar=ch$ystar[s-1,,drop=FALSE])
        snapshot$Vcomps<<-get('Vcomps',env);snapshot$Z<<-get('Z',env)
      }
      if(length(fidx))snapshot$chain_index<<-get('i',frames[[tail(fidx,1)]])
      stopifnot(identical(rng,get('.Random.seed',envir=.GlobalEnv)))
    },error=function(ex)inspection_error<<-conditionMessage(ex))
  }),error=identity)
saveRDS(list(snapshot=snapshot,inspection_error=inspection_error,
  outcome=if(inherits(result,'error'))conditionMessage(result) else 'completed'),
  file.path(out,'failure_case/error_time_trace.rds'))
stopifnot(inherits(result,'error'),is.null(inspection_error),!is.null(snapshot))
capture.output(str(snapshot[c('message','conditionCall','iteration','chain_index','previous_state')],
  max.level=2),file=file.path(out,'failure_case/trace_summary.txt'))
writeLines(unlist(snapshot$calls),file.path(out,'failure_case/call_stack.txt'))
matrix_checks<-function(A,label,provenance){
  finite<-all(is.finite(A));ev<-if(finite)eigen((A+t(A))/2,symmetric=TRUE,only.values=TRUE)$values else NA_real_
  ch<-tryCatch({chol(A);TRUE},error=function(e)FALSE)
  determinant_info<-if(finite)determinant(A,logarithm=TRUE) else list(sign=NA_real_,modulus=NA_real_)
  data.frame(matrix=label,provenance=provenance,dimension=nrow(A),nonfinite=sum(!is.finite(A)),
    max_asymmetry=if(finite)max(abs(A-t(A))) else NA_real_,
    symmetrized_min_eigenvalue=min(ev),symmetrized_max_eigenvalue=max(ev),cholesky_success=ch,
    determinant=if(finite)det(A) else NA_real_,determinant_sign=determinant_info$sign,
    log_abs_determinant=as.numeric(determinant_info$modulus))
}
H<-snapshot$rtmvnorm_arguments$H
checks<-matrix_checks(H,'H','actual rtmvnorm argument from failing call')
saveRDS(H,file.path(out,'failure_case/failing_H.rds'))
if(!is.null(snapshot$previous_state)){
  # K and V reconstructed AFTER the failure, using recorded values and evaluated implementation.
  K<-exp(-getFromNamespace('makeKpart','bkmr')(as.numeric(snapshot$previous_state$r),snapshot$Z))
  V<-diag(nrow(K))+as.numeric(snapshot$previous_state$lambda)*K
  checks<-rbind(checks,matrix_checks(K,'K','reconstructed from recorded r and Z after failure'),
                matrix_checks(V,'V','reconstructed I + lambda K after failure'))
  saveRDS(list(K=K,V=V),file.path(out,'failure_case/reconstructed_K_V.rds'))
}
write.csv(checks,file.path(out,'failure_case/matrix_checks.csv'),row.names=FALSE)
capture.output(print(tmvtnorm::rtmvnorm),file=file.path(out,'audit/tmvtnorm_rtmvnorm.R.txt'))
capture.output(print(getFromNamespace('checkSymmetricPositiveDefinite','tmvtnorm')),
  file=file.path(out,'audit/tmvtnorm_checkSymmetricPositiveDefinite.R.txt'))
print(checks)
