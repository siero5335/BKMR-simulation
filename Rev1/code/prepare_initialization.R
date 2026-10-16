# Run from an isolated prepared workspace, after the 12 baseline refits.
source('code/validate_inputs.R')
out<-'initialization_extension_20260914'
args<-commandArgs(TRUE)
ids<-if(length(args))args else names(specs)
stopifnot(all(ids %in% names(specs)))
for(id in ids){
 s<-specs[[id]]
 cp<-readRDS(sprintf('revision_20260911/checkpoints/task_%03d.rds',as.integer(id)))
 stopifnot(cp$current_status=='ok',isTRUE(cp$identity_matches),identical(cp$starting_values,s$initial_sets$original))
 stopifnot(identical(tools::md5sum(names(cp$fingerprint$packages)),cp$fingerprint$packages))
 inp<-list(task=s$task,data=s$data,initial_sets=s$initial_sets,reference_draws=cp$draws,
           baseline_diagnostics=s$baseline_diagnostics,fingerprint=cp$fingerprint)
 inp$fingerprint$RNGkind<-s$RNGkind
 p<-file.path(out,'inputs',sprintf('task_%03d.rds',as.integer(id)))
 if(file.exists(p))stopifnot(identical(readRDS(p),inp)) else saveRDS(inp,p)
}
cat('Created extension inputs with regenerated baseline reference draws. Historical trajectories were not distributed.\n')
