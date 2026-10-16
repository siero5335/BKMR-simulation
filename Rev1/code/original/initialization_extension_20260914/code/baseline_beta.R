options(stringsAsFactors=FALSE)
out<-'initialization_extension_20260914'
reg<-read.csv('h_false_allpoints_20260912/registry.csv');by_chain<-list();by_task<-list();j<-0L;k<-0L
for(i in seq_len(nrow(reg))){
 x<-readRDS(reg$checkpoint[i]);if(x$status!='completed')next
 m<-x$draws[['beta[1]']];z<-x$diagnostics[x$diagnostics$parameter=='beta[1]',]
 cm<-colMeans(m);cq<-apply(m,2,quantile,probs=c(.025,.975))
 for(ch in 1:4){j<-j+1L;by_chain[[j]]<-data.frame(task_index=x$task$task_index,n=x$task$n,chain=ch,
  sampler_seed=as.integer(strsplit(x$task$chain_seeds,';',fixed=TRUE)[[1]])[ch],
  beta_start=x$starting_values[[ch]]$beta,lambda_start=x$starting_values[[ch]]$lambda,
  beta_mean=cm[ch],beta_q025=cq[1,ch],beta_q975=cq[2,ch],beta_error=cm[ch]-.5)}
 k<-k+1L;by_task[[k]]<-data.frame(task_index=x$task$task_index,n=x$task$n,primary_pass=x$primary_pass,
  beta_rhat=z$rhat,beta_bulk_ess=z$bulk_ess,beta_tail_ess=z$tail_ess,beta_mean=z$mean,beta_sd=z$sd,
  beta_q025=z$q025,beta_q975=z$q975,chain_mean_min=min(cm),chain_mean_max=max(cm),
  chain_mean_range=max(cm)-min(cm),chain_mean_range_over_pooled_sd=(max(cm)-min(cm))/z$sd)
}
write.csv(do.call(rbind,by_chain),file.path(out,'baseline_beta_by_chain.csv'),row.names=FALSE)
write.csv(do.call(rbind,by_task),file.path(out,'baseline_beta_by_task.csv'),row.names=FALSE)
cat('Baseline: 199 datasets and 796 chain-specific beta summaries.\n')
