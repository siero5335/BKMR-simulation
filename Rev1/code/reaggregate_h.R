# Base-R reaggregation of all reported point-level diagnostics; no raw-draw recalculation.
script<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))][1])
root<-normalizePath(file.path(dirname(script),'..'))
h<-read.csv(file.path(root,'h_observed_points/h_by_task_point.csv'))
core<-read.csv(file.path(root,'reference/h_core_diagnostics_199.csv'))
task<-read.csv(file.path(root,'h_observed_points/by_task.csv'))
meets<-function(d)is.finite(d$rhat)&is.finite(d$bulk_ess)&is.finite(d$tail_ess)&d$rhat<=1.01&d$bulk_ess>=400&d$tail_ess>=400
h$recomputed_pass<-meets(h)
stopifnot(nrow(h)==29680,!anyDuplicated(paste(h$task_index,h$observed_index)),nrow(core)==1194)
rows<-lapply(split(h,h$task_index),function(g){
 id<-g$task_index[1];n<-g$n[1]
 stopifnot(identical(sort(g$observed_index),seq_len(n)))
 c<-core[core$task_index==id,];stopifnot(nrow(c)==6)
 primary<-all(meets(c));original<-task[task$task_index==id,]
 stopifnot(nrow(original)==1,primary==(tolower(original$primary_pass)=='true'))
 five<-g$observed_index %in% c(1,n/4,n/2,3*n/4,n)
 data.frame(task_index=id,n,primary_pass=primary,series=n,series_meeting=sum(g$recomputed_pass),
  all_points=all(g$recomputed_pass),nested_five=all(g$recomputed_pass[five]))
})
r<-do.call(rbind,rows)
s<-do.call(rbind,lapply(split(r,r$n),function(g)data.frame(n=g$n[1],datasets=nrow(g),
 h_series=sum(g$series),series_meeting=sum(g$series_meeting),datasets_all_h_meet=sum(g$all_points),
 five_meet_other_flagged=sum(g$nested_five & !g$all_points))))
expected<-read.csv(file.path(root,'h_observed_points/h_summary.csv'));expected<-expected[expected$scope=='all observed points',]
for(i in seq_len(nrow(s))){e<-expected[expected$n==s$n[i],];for(k in c('datasets','h_series','series_meeting','datasets_all_h_meet'))stopifnot(s[[k]][i]==e[[k]])}
stopifnot(nrow(r)==199,sum(r$primary_pass)==22,sum(r$all_points & !r$primary_pass)==84,
 sum(r$all_points & r$primary_pass)==22,sum(r$nested_five)==119,sum(r$nested_five & !r$all_points)==13)
args<-commandArgs(TRUE)
dest<-if(length(args))args[1] else file.path(root,'_work/h_checks')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write.csv(s,file.path(dest,'independent_all_point_aggregation.csv'),row.names=FALSE)
cat('PASS: all 29,680 diagnostic records reaggregated; primary mapping and five/all comparison agree.\n');print(s)
