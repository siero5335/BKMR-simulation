out<-'initialization_extension_20260914'
f<-read.csv(file.path(out,'fit_summaries.csv'));d<-read.csv(file.path(out,'parameter_diagnostics.csv'))
cases<-read.csv(file.path(out,'cases.csv'))
cols<-c(original='#2455A4',alternative_1='#C05B19',alternative_2='#22845B')
budgets<-c(10000,20000,40000,100000);xx<-log10(budgets)
title_for<-function(row)paste0('n=',row$n,'; task ',row$task_index,' (',
 if(row$selection_role=='previous_anchor')'anchor' else if(row$selection_role=='lowest_id_primary_pass')'pass control' else if(grepl('middle',row$selection_role))'median non-pass' else 'high beta Rhat',')')
for(kind in c('beta_intervals','maximum_rhat','minimum_ess')){
 png(file.path(out,'figures',paste0(kind,'.png')),width=2100,height=2300,res=190)
 par(mfrow=c(4,3),mar=c(3.2,4.1,2.6,.5),oma=c(3,0,4,0))
 for(i in seq_len(nrow(cases))){
  id<-cases$task_index[i];ff<-f[f$task_index==id & f$window=='fixed_5000',]
  linesets<-list()
  for(set in names(cols)){
   q<-ff[ff$initial_set==set,];q<-q[match(budgets,q$budget),]
   if(kind=='beta_intervals')v<-q$beta_mean else {
    v<-vapply(budgets,function(b){
     z<-d[d$task_index==id & d$initial_set==set & d$window=='fixed_5000' & d$budget==b,]
     vals<-if(kind=='maximum_rhat')z$rhat else c(z$bulk_ess,z$tail_ess)
     if(nrow(z)!=6 || any(!is.finite(vals)))NA_real_ else if(kind=='maximum_rhat')max(vals) else min(vals)
    },numeric(1))
   }
   linesets[[set]]<-list(q=q,v=v)
  }
  vals<-unlist(lapply(linesets,function(z)if(kind=='beta_intervals')c(z$q$beta_q025,z$q$beta_q975,.5) else z$v))
  vals<-vals[is.finite(vals)]
  if(kind=='maximum_rhat')vals<-c(vals,.99,1.01)
  if(kind=='minimum_ess')vals<-c(vals,400)
  yl<-if(length(vals))range(vals) else c(0,1)
  if(diff(yl)==0)yl<-yl+c(-.01,.01)
  if(kind=='minimum_ess')yl[1]<-max(.01,yl[1]*.8)
  plot(xx,rep(NA_real_,4),type='n',xlim=range(xx)+c(-.05,.05),ylim=yl,
   log=if(kind %in% c('minimum_ess','maximum_rhat'))'y' else '',xaxt='n',
   xlab='',ylab=if(kind=='beta_intervals')'Beta mean and 95% interval' else if(kind=='maximum_rhat')'Maximum Rhat (six parameters)' else 'Minimum bulk/tail ESS (six parameters)',
   main=title_for(cases[i,]),cex.main=.85)
  axis(1,at=xx,labels=c('10k','20k','40k','100k'),cex.axis=.8)
  abline(h=if(kind=='beta_intervals').5 else if(kind=='maximum_rhat')1.01 else 400,lty=2,col='grey45')
  for(j in seq_along(cols)){
   set<-names(cols)[j];z<-linesets[[set]];xp<-xx+(j-2)*.012
   if(kind=='beta_intervals')segments(xp,z$q$beta_q025,xp,z$q$beta_q975,col=adjustcolor(cols[j],alpha.f=.5))
   lines(xp,z$v,col=cols[j],lwd=1.4)
   points(xp,z$v,col=cols[j],pch=ifelse(is.na(z$q$primary_pass),4,ifelse(z$q$primary_pass,16,1)),cex=.8)
  }
 }
 mtext(switch(kind,beta_intervals='Initialization and iteration sensitivity: beta',maximum_rhat='Initialization and iteration sensitivity: Rhat',minimum_ess='Initialization and iteration sensitivity: effective sample size'),outer=TRUE,side=3,line=2,cex=1.2)
 mtext('Blue: original starts; orange/green: alternative starts. Filled points: six-parameter primary pass. Warmup fixed at 5,000.',outer=TRUE,side=1,line=1,cex=.8)
 mtext('Missing points indicate unavailable or non-finite summaries. Cases were selected using baseline diagnostics.',outer=TRUE,side=1,line=2,cex=.8)
 dev.off()
}
bc<-read.csv(file.path(out,'beta_by_chain.csv'))
png(file.path(out,'figures','beta_chain_means.png'),width=2100,height=2300,res=190)
par(mfrow=c(4,3),mar=c(3.2,4.1,2.6,.5),oma=c(3,0,4,0))
for(i in seq_len(nrow(cases))){
 id<-cases$task_index[i]
 z<-bc[bc$task_index==id & bc$window=='fixed_5000' & bc$budget %in% c(10000,100000),]
 yl<-range(c(z$beta_mean,.5),finite=TRUE)
 if(diff(yl)==0)yl<-yl+c(-.01,.01)
 plot(NA,xlim=c(.5,12.5),ylim=yl,xaxt='n',xlab='',ylab='Chain-specific beta mean',
  main=title_for(cases[i,]),cex.main=.85)
 abline(h=.5,lty=2,col='grey45');abline(v=c(4.5,8.5),col='grey85')
 axis(1,at=c(2.5,6.5,10.5),labels=c('Original','Alternative 1','Alternative 2'),cex.axis=.8)
 for(j in seq_along(cols))for(ch in 1:4){
  a<-z[z$initial_set==names(cols)[j] & z$chain==ch & z$budget==10000,]
  b<-z[z$initial_set==names(cols)[j] & z$chain==ch & z$budget==100000,]
  xp<-(j-1)*4+ch
  if(nrow(a)==1 && nrow(b)==1)segments(xp,a$beta_mean,xp,b$beta_mean,col=cols[j],lwd=1.3)
  if(nrow(a)==1)points(xp,a$beta_mean,col=cols[j],pch=1,cex=.85)
  if(nrow(b)==1)points(xp,b$beta_mean,col=cols[j],pch=16,cex=.8)
 }
}
mtext('Initial states and chain-specific beta estimates',outer=TRUE,side=3,line=2,cex=1.2)
mtext('Within each group, chains 1 to 4 appear from left to right. Open circles: 10k iterations; filled circles: 100k.',outer=TRUE,side=1,line=1,cex=.8)
mtext('Warmup fixed at 5,000. Dashed line: generating beta = 0.5. Groups use paired sampler seeds; do not pool them.',outer=TRUE,side=1,line=2,cex=.8)
dev.off()
