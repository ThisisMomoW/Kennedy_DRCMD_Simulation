kennedy_score <- function(dat,kappa,lambda1,gamma1,beta0,beta1) {
 n<-nrow(dat);stopifnot(length(kappa)==n,length(lambda1)==n,length(gamma1)==n,
  length(beta0)==n,length(beta1)==n,all(kappa>0),all(gamma1>0 & gamma1<1))
 observed<-which(dat$R==1);phi<-matrix(NA_real_,n,2,dimnames=list(NULL,c('phi0','phi1')))
 for(a in 0:1) {
  gamma<-if(a==1) gamma1 else 1-gamma1;lambda<-if(a==1) lambda1 else 1-lambda1
  m<-(if(a==1) beta1 else beta0)/gamma
  residual<-numeric(n) # Explicitly zero for R=0: never multiply NA by zero.
  residual[observed]<-((dat$A_obs[observed]==a)-lambda[observed])/kappa[observed]
  phi[,a+1]<-m+(dat$Y-m)/gamma*(lambda+residual)
 }
 phi
}
kennedy_finish <- function(phi) {
 delta<-phi[,2]-phi[,1];estimate<-mean(delta);se<-sd(delta)/sqrt(nrow(phi))
 list(estimate=estimate,se=se,lower=estimate-1.96*se,upper=estimate+1.96*se,
  risk0=mean(phi[,1]),risk1=mean(phi[,2]),point_success=is.finite(estimate),
  inference_success=all(is.finite(c(estimate,se))) && se>=0)
}
kennedy_fold <- function(train,evaluate,train_ids,test_ids) {
 stopifnot(length(intersect(train_ids,test_ids))==0,nrow(train)==length(train_ids),nrow(evaluate)==length(test_ids))
 training<-train[,c('X1','X2','Y','R','A_obs'),drop=FALSE]
 eval_features<-evaluate[,c('X1','X2','Y'),drop=FALSE]
 stages<-list();models<-list()
 sl<-function(y,x,family,ids,stage) {
  stopifnot(length(y)==nrow(x),length(ids)==length(y),!anyNA(y),all(ids %in% train_ids),!any(ids %in% test_ids))
  fit<-SuperLearner::SuperLearner(Y=as.numeric(y),X=x,family=family,SL.library='SL.gam',
    obsWeights=rep(1,length(y)),cvControl=list(V=5),env=asNamespace('SuperLearner'))
  # These fields contain training-only results, so signatures can verify fold isolation.
  stages[[stage]]<<-list(stage=stage,family=family$family,library='SL.gam',fit_ids=ids,
   rows=length(y),cv_requested=5L,cv_actual=length(fit$validRows),cv_validation_rows=fit$validRows,
   coef=fit$coef,cv_risk=fit$cvRisk,
   training_signature=object_md5(list(y=y,x=x,coef=fit$coef,cvRisk=fit$cvRisk,Z=fit$Z,validRows=fit$validRows)))
  models[[stage]]<<-fit
  fit
 }
 pred<-function(fit,x) as.numeric(predict(fit,newdata=x)$pred)
 xz<-training[,c('X1','X2','Y'),drop=FALSE];xx<-training[,c('X1','X2'),drop=FALSE]
 ez<-eval_features;ex<-eval_features[,c('X1','X2'),drop=FALSE];obs<-which(training$R==1)
 # SL.gam uses response symbol Y internally; rename this predictor only, preserving values.
 names(xz)[names(xz)=='Y']<-'outcome_Y';names(ez)[names(ez)=='Y']<-'outcome_Y'
 kfit<-sl(training$R,xz,binomial(),train_ids,'kappa')
 lfit<-sl(training$A_obs[obs],xz[obs,,drop=FALSE],binomial(),train_ids[obs],'lambda')
 lambda_train<-pred(lfit,xz);lambda_test<-pred(lfit,ez)
 gfit<-sl(lambda_train,xx,gaussian(),train_ids,'gamma1')
 b1fit<-sl(training$Y*lambda_train,xx,gaussian(),train_ids,'beta1')
 b0fit<-sl(training$Y*(1-lambda_train),xx,gaussian(),train_ids,'beta0')
 raw<-data.frame(row=test_ids,kappa_raw=pred(kfit,ez),lambda1=lambda_test,gamma1_raw=pred(gfit,ex),
  beta0=pred(b0fit,ex),beta1=pred(b1fit,ex))
 raw$kappa<-pmin(.99,pmax(.01,raw$kappa_raw));raw$gamma1<-pmin(.99,pmax(.01,raw$gamma1_raw));raw$gamma0<-1-raw$gamma1
 raw$m0<-raw$beta0/raw$gamma0;raw$m1<-raw$beta1/raw$gamma1
 phi<-kennedy_score(evaluate,raw$kappa,raw$lambda1,raw$gamma1,raw$beta0,raw$beta1)
 list(phi=phi,predictions=raw,fit_audit=stages,
  generated_training=data.frame(row=train_ids,lambda1=lambda_train,gamma_response=lambda_train,
   beta1_response=training$Y*lambda_train,beta0_response=training$Y*(1-lambda_train)),
  truncation=data.frame(n=length(test_ids),kappa=sum(raw$kappa!=raw$kappa_raw,na.rm=TRUE),
   gamma1=sum(raw$gamma1!=raw$gamma1_raw,na.rm=TRUE)),
  train_ids=train_ids,test_ids=test_ids)
}
kennedy_fit <- function(dat,folds) {
 # Select only observed variables even if the caller includes any evaluation-only data.
 dat<-dat[,c('X1','X2','Y','R','A_obs'),drop=FALSE]
 stopifnot(all(is.na(dat$A_obs)==(dat$R==0)))
 verify_folds(folds,nrow(dat));pieces<-list();phi<-matrix(NA_real_,nrow(dat),2)
 for(k in seq_along(folds)) {
  f<-folds[[k]];piece<-kennedy_fold(dat[f$train,,drop=FALSE],dat[f$test,,drop=FALSE],f$train,f$test)
  pieces[[k]]<-piece;phi[f$test,]<-piece$phi
 }
 list(result=kennedy_finish(phi),scores=phi,fold_details=pieces)
}
verify_folds <- function(folds,n) {
 stopifnot(length(folds)==5L,identical(sort(as.integer(unlist(lapply(folds,`[[`,'test')))),seq_len(n)))
 for(f in folds) stopifnot(length(intersect(f$train,f$test))==0,
  identical(sort(as.integer(c(f$train,f$test))),seq_len(n)))
 invisible(TRUE)
}
