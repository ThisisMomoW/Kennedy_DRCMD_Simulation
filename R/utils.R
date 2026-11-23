atomic_save <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path, '.', Sys.getpid(), '.tmp')
  saveRDS(x, tmp, version=3)
  stopifnot(file.rename(tmp,path))
}
object_md5 <- function(x) {
  f<-tempfile();on.exit(unlink(f));saveRDS(x,f,version=3);unname(tools::md5sum(f))
}
set_stream <- function(s) assign('.Random.seed',s,.GlobalEnv)
close_numeric <- function(x,y,tol=1e-10) {
  stopifnot(length(x)==length(y),identical(is.na(x),is.na(y)))
  ok<-!is.na(x);stopifnot(all(abs(x[ok]-y[ok])<=tol*(1+abs(y[ok]))));invisible(TRUE)
}
reset_audit <- function() {
  e<-get('.gam_audit',.GlobalEnv)
  e$counts<-setNames(integer(7),c('drcmd','est_kappa','est_g','est_m_a','est_varphi','est_varphi_eem','est_psi_tml'))
  e$stage<-'';e$folds<-NULL;e$sl<-e$cv<-e$clipping<-list()
}
read_audit <- function() as.list(get('.gam_audit',.GlobalEnv))
scenario_paths<-c('scenarios/01_observation_heterogeneity','scenarios/02_outcome_slope','scenarios/03_nonlinear_outcome')
method_ids<-c('original','direct_m_gam','drcmd')
method_labels<-c(original='Kennedy ratio GAM',direct_m_gam='Kennedy direct GAM',drcmd='DR-CMD')
software_info<-function() {
  pkgs<-c('drcmd','SuperLearner','gam','nnls','ggplot2','patchwork','rmarkdown','knitr')
  do.call(rbind,lapply(pkgs,function(p){d<-packageDescription(p)
    data.frame(package=p,version=d$Version,commit=if(is.null(d$RemoteSha))''else d$RemoteSha,path=find.package(p))}))
}
