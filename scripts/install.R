# Optional setup in a new checkout. No package is upgraded during simulation.
dir.create('library',showWarnings=FALSE)
.libPaths(unique(c(normalizePath('library'),.libPaths())))
options(repos=c(CRAN='https://cloud.r-project.org'))
if(!requireNamespace('remotes',quietly=TRUE))install.packages('remotes',lib='library')
versions<-read.csv('software/versions.csv')
for(i in which(versions$package!='drcmd')){
 p<-versions$package[i];v<-versions$version[i]
 if(!requireNamespace(p,quietly=TRUE)||packageVersion(p)!=package_version(v))
  remotes::install_version(p,version=v,lib='library',upgrade='never',dependencies=NA)
}
# Official package source and commit metadata are included in software/drcmd.
install.packages('software/drcmd',repos=NULL,type='source',lib='library')
stopifnot(packageDescription('drcmd',lib.loc='library')$RemoteSha=='76dc26abdb921025a2710b2826e96efe9b9df800')
