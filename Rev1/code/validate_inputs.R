# Generate inputs/starts only; does not run MCMC.
source('Main/Code_Data/code/_pure_bkmr_core.R')
specs<-dget('reference/initialization_specs.R')
tasks<-read.csv('Main/Code_Data/results/task_table.csv')
base<-read.csv('Main/Code_Data/results/all_task_results.csv')
reg<-read.csv('initialization_extension_20260914/registry.csv')
for(id in names(specs)){
 s<-specs[[id]];t<-tasks[tasks$task_index==as.integer(id),]
 do.call(RNGkind,as.list(s$RNGkind))
 dat<-make_bkmr_dataset(t$dataset_seed,t$n,t$hfun,t$beta_true)
 stopifnot(identical(dat,s$data))
 seed<-derive_initial_value_seed(t$initial_value_seed,t$scenario,t$dataset_seed,t$n,t$hfun,t$beta_true,t$iter,t$warmup,t$varsel,t$chain_seeds)
 stopifnot(seed==s$initial_value_task_seed,dat$label_sum==base$label_sum[base$task_index==as.integer(id)])
 for(set in c('original','alternative_1','alternative_2')){
  g<-reg[reg$task_index==as.integer(id)&reg$initial_set==set,]
  stopifnot(length(unique(g$initial_seed))==1)
  starts<-make_bkmr_binomial_starting_values(dat$y,dat$Z,dat$X,4,FALSE,g$initial_seed[1])
  stopifnot(identical(starts,s$initial_sets[[set]]))
 }
}
cat('PASS: all 12 simulated datasets and 144 initial states exactly match exported historical inputs; no sampling.\n')
