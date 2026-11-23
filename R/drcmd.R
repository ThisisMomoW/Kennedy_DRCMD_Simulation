# Native package fitting and the arm-difference risk-difference output.
drcmd_point <- function(dat) {
 stopifnot(identical(is.na(dat$A_obs),dat$R==0));reset_audit()
 fit<-capture_fit(drcmd::drcmd(Y=dat$Y,A=dat$A_obs,X=dat[,c('X1','X2'),drop=FALSE],
  default_learners='SL.gam',eem_ind=FALSE,tml=FALSE,k=5,cv_folds=5,cutoff=.01,parallel=FALSE))
 v<-fit$value; est<-se<-NA_real_; risks<-c(NA_real_,NA_real_)
 if(!is.null(v)) {
  stopifnot(all(is.na(v$params$Rprobs)),identical(v$U,'A'))
  est<-as.numeric(v$results$estimates$psi_hat_ate);se<-as.numeric(v$results$ses$psi_hat_ate)
  risks<-c(v$results$estimates$psi_0_hat,v$results$estimates$psi_1_hat)
 }
 list(estimate=est,se=se,lower=est-1.96*se,upper=est+1.96*se,risks=risks,fit=fit,audit=read_audit(),
  inference_success=all(is.finite(c(est,se))) && se>=0)
}

