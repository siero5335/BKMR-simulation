args<-commandArgs(TRUE)
stopifnot(length(args)==1)
versions<-read.csv(args[1])
for(i in seq_len(nrow(versions))){
 p<-versions$package[i]
 if(!requireNamespace(p,quietly=TRUE))stop('Missing package: ',p)
 if(as.character(packageVersion(p))!=versions$version[i])stop('Version mismatch for ',p,': expected ',versions$version[i])
}
stopifnot(identical(packageDescription('bkmr')$RemoteSha,'45413e338a316362d629f53bd2a917c4bf485c1e'))
print(sessionInfo());print(RNGkind())
cat('Required package versions and bkmr commit match. This does not establish binary/BLAS identity.\n')
