out <- 'initialization_extension_20260914'
reg <- read.csv(file.path(out,'registry.csv'),stringsAsFactors=FALSE)
rows <- lapply(seq_len(nrow(reg)),function(i){
 r <- reg[i,]
 inp <- readRDS(file.path(out,'inputs',sprintf('task_%03d.rds',r$task_index)))
 v <- inp$initial_sets[[r$initial_set]][[r$chain]]
 data.frame(job_id=r$job_id,task_index=r$task_index,n=r$n,initial_set=r$initial_set,
  chain=r$chain,sampler_seed=r$sampler_seed,beta_start=v$beta,lambda_start=v$lambda,
  r_z1_start=v$r[1],r_z2_start=v$r[2],r_z3_start=v$r[3],r_z4_start=v$r[4],
  h_hat_start_min=min(v$h.hat),h_hat_start_max=max(v$h.hat),
  ystar_start_min=min(v$ystar),ystar_start_max=max(v$ystar))
})
write.csv(do.call(rbind,rows),file.path(out,'initial_values.csv'),row.names=FALSE)
