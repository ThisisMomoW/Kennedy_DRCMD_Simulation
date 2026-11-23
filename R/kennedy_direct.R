# Arm-specific outcome GAMs reuse the fold-specific exposure-imputation weights.
# Use the installed SL.gam wrapper, not a reimplementation or mgcv substitute.
dmg_start_audit<-function() {
 assign('.dmg_audit',new.env(parent=emptyenv()),.GlobalEnv)
 suppressMessages(trace('gam.fit',tracer=quote({
  ee<-get('.dmg_audit',.GlobalEnv);ee$input_weights<-as.numeric(weights)
  ee$family<-family$family;ee$link<-family$link
 }),exit=quote({
  ee<-get('.dmg_audit',.GlobalEnv)
  ee$termination<-list(iterations=iter,maxit=maxit,backfit_iterations=bf.iter,
   old_deviance=old.dev,new_deviance=new.dev,epsilon=epsilon,
   needs_another_iteration=one.more,
   criterion_met=is.finite(new.dev)&&is.finite(old.dev)&&!isTRUE(one.more),
   deviance_relative_change=abs(old.dev-new.dev)/(old.dev+.1))
 }),where=asNamespace('gam'),print=FALSE))
}
dmg_fit_fold<-function(training,evaluation,lambda_train) {
 tr<-training[,c('X1','X2','Y'),drop=FALSE];xx<-tr[,c('X1','X2'),drop=FALSE]
 ex<-evaluation[,c('X1','X2'),drop=FALSE]
 stopifnot(length(lambda_train)==nrow(tr),all(is.finite(lambda_train)),all(lambda_train>=0&lambda_train<=1),
  all(vapply(xx,function(x)length(unique(x))>4,TRUE)))
 fits<-list();pred<-vector('list',2)
 for(a in 0:1) {
  w<-if(a==1)lambda_train else 1-lambda_train
  ee<-get('.dmg_audit',.GlobalEnv);rm(list=ls(ee),envir=ee)
  cap<-capture_fit({
   z<-SuperLearner::SL.gam(Y=tr$Y,X=xx,newX=ex,family=quasibinomial(link='logit'),obsWeights=w,deg.gam=2,cts.num=4)
   m<-z$fit$object;pp<-as.numeric(z$pred)
   stopifnot(all(is.finite(pp)),length(pp)==nrow(ex),length(m$prior.weights)==length(w),
    identical(as.numeric(m$prior.weights),as.numeric(w)),identical(ee$input_weights,as.numeric(w)),
    m$family$family=='quasibinomial',m$family$link=='logit',ee$family=='quasibinomial',ee$link=='logit',
    identical(attr(terms(m),'term.labels'),c('s(X1, 2)','s(X2, 2)')))
   list(prediction=pp,coefficients=coef(m),formula=paste(deparse(formula(m)),collapse=' '),
    prior_weights=as.numeric(m$prior.weights),model_rows=nrow(m$model),control=m$control,
    nonlinear_df=m$nl.df,df_residual=m$df.residual,rank=m$rank,
    all_iter_fields=m[which(names(m)=='iter')],backfit_iterations=m$bf.iter,deviance=m$deviance,
    convergence=ee$termination)
  })
  # Preserve failure-stage instrumentation even if fitting/prediction fails.
  cap$gam_fit_audit<-as.list(ee);cap$weights<-w;cap$arm<-a
  fits[[a+1L]]<-cap;pred[a+1L]<-list(cap$value$prediction)
 }
 list(fits=fits,m0=pred[[1]],m1=pred[[2]],success=all(vapply(fits,function(x)!is.null(x$value),TRUE)))
}
