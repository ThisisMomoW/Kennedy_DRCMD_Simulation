# Saved-fit diagnostics and controlled, seed-matched replay of failed Kennedy fits.
source('R/load.R');m<-check_study();r<-read.csv('results/replicates.csv');s<-read.csv('results/summary.csv')
stopifnot(nrow(r)==9600)
dir.create('results/diagnostics',showWarnings=FALSE)
index_for<-function(key)match(key,m$jobs$key)
path_for<-function(key,part)checkpoint_path(m$jobs[index_for(key),],part)
# For each cell and method, inspect its largest squared error without excluding it.
worst<-do.call(rbind,lapply(split(subset(r,point_valid),interaction(subset(r,point_valid)$scenario,subset(r,point_valid)$method,drop=TRUE)),function(z)z[which.max((z$estimate-z$truth)^2),]))
out<-obs_rows<-weight_rows<-list();oi<-wi<-0L
for(i in seq_len(nrow(worst))){
 z<-worst[i,];dat<-readRDS(path_for(z$key,'data'))$observed;folds<-readRDS(path_for(z$key,'folds'))
 if(z$method=='drcmd'){
  x<-readRDS(path_for(z$key,'drcmd'));scores<-rep(NA_real_,nrow(dat))
  for(k in seq_along(folds)){
   v<-x$audit$fold_results[[k]];cc<-x$audit$centering[[k]]
   # psi_ate_ic belongs to the unused direct-contrast fit. Reconstruct the
   # reported arm difference from the two arm corrections and training constants.
   scores[v$idx]<-v$ics$psi_1_ic-v$ics$psi_0_ic+cc$plugin1-cc$plugin0
  }
  close_numeric(mean(scores),z$estimate)
  for(v in x$audit$sl){w<-v$weights
   wi<-wi+1L;weight_rows[[wi]]<-data.frame(key=z$key,method=z$method,fold=v$fold,stage=v$stage,role=v$role,
     n=length(w),positive=sum(w>0),maximum=max(w),p99=unname(quantile(w,.99)),ess=sum(w)^2/sum(w^2),family=v$family)
  }
  pos<-head(order(abs(scores),decreasing=TRUE),10)
  obs<-data.frame(row=pos,score=scores[pos]);ratio_min<-ratio_max<-gamma_clips<-NA_real_
 }else{
  kr<-readRDS(path_for(z$key,'original'))$value
  pr<-do.call(rbind,lapply(seq_along(folds),function(k)cbind(fold=k,kr$fold_details[[k]]$predictions)))
  pr<-pr[order(pr$row),];stopifnot(identical(as.integer(pr$row),seq_len(nrow(dat))))
  if(z$method=='direct_m_gam'){
   kd<-readRDS(path_for(z$key,'direct_m_gam'));pr$m0<-pr$m1<-NA_real_
   for(k in seq_along(folds)){
    pr$m0[folds[[k]]$test]<-kd$fold_details[[k]]$m0;pr$m1[folds[[k]]$test]<-kd$fold_details[[k]]$m1
   }
   expected<-kd$scores
  }else expected<-kr$scores
  corrections<-matrix(0,nrow(dat),2);phi<-corrections;observed<-which(dat$R==1)
  for(a in 0:1){ga<-if(a)pr$gamma1 else pr$gamma0;la<-if(a)pr$lambda1 else 1-pr$lambda1;ma<-if(a)pr$m1 else pr$m0
   residual<-numeric(nrow(dat));residual[observed]<-((dat$A_obs[observed]==a)-la[observed])/pr$kappa[observed]
   corrections[,a+1]<-(dat$Y-ma)/ga*(la+residual);phi[,a+1]<-ma+corrections[,a+1]
  }
  stopifnot(max(abs(phi-expected))<1e-9*(1+max(abs(expected))))
  scores<-phi[,2]-phi[,1];close_numeric(mean(scores),z$estimate)
  pos<-head(order(abs(scores),decreasing=TRUE),10)
  obs<-data.frame(pr[pos,],Y=dat$Y[pos],R=dat$R[pos],A_obs=dat$A_obs[pos],
                  correction0=corrections[pos,1],correction1=corrections[pos,2],score=scores[pos])
  ratio_min<-min(pr$m0,pr$m1);ratio_max<-max(pr$m0,pr$m1);gamma_clips<-sum(pr$gamma1_raw!=pr$gamma1)
 }
 # Individual-score contributions are signed; squared-error shares are setting-level summaries.
 out[[i]]<-data.frame(z[,c('key','scenario','group','n','level','method','estimate','truth')],
   max_abs_individual_score=max(abs(scores)),top10_signed_contribution=sum(scores[pos])/nrow(dat),
   remainder_signed_contribution=sum(scores[-pos])/nrow(dat),mean_prediction_min=ratio_min,mean_prediction_max=ratio_max,gamma_truncations=gamma_clips)
 obs$key<-z$key;obs$method<-z$method
 write.csv(obs,file.path('results/diagnostics',paste0(z$key,'_',z$method,'_largest_scores.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,out),'results/diagnostics/worst_fit_components.csv',row.names=FALSE)
write.csv(do.call(rbind,weight_rows),'results/diagnostics/drcmd_weights.csv',row.names=FALSE)
# A replay is diagnostic only: it never replaces an attempt's saved result.
failed<-subset(r,method=='original'&!point_valid);replays<-list()
for(j in seq_len(nrow(failed))){
 z<-failed[j,];dat<-readRDS(path_for(z$key,'data'))$observed;folds<-readRDS(path_for(z$key,'folds'))
 seed<-m$seeds[[index_for(z$key)]]$kennedy
 assign('.failure_trace',new.env(parent=emptyenv()),.GlobalEnv);e<-get('.failure_trace',.GlobalEnv)
 e$stages<-list();e$errors<-list();e$fold<-0L;e$stage<-'';e$calls<-list();e$events<-list()
 suppressMessages(trace('kennedy_fold',tracer=quote({e<-get('.failure_trace',.GlobalEnv);e$fold<-e$fold+1L;e$stage_count<-0L;e$train<-train;e$test_ids<-test_ids}),print=FALSE))
 suppressMessages(trace('SuperLearner',tracer=quote({
   e<-get('.failure_trace',.GlobalEnv);e$stage_count<-e$stage_count+1L;e$stage<-c('kappa','lambda','gamma1','beta1','beta0')[e$stage_count];e$gam_count<-0L
   e$stages[[length(e$stages)+1L]]<-list(fold=e$fold,stage=e$stage,n=length(Y),response_range=range(Y),family=family$family,
     zeros=sum(Y==0),ones=sum(Y==1),observed0=sum(e$train$R==1&e$train$A_obs==0,na.rm=TRUE),observed1=sum(e$train$R==1&e$train$A_obs==1,na.rm=TRUE))
 }),where=asNamespace('SuperLearner'),print=FALSE))
 suppressMessages(trace('SL.gam',tracer=quote({
   e<-get('.failure_trace',.GlobalEnv);e$gam_count<-e$gam_count+1L;e$gam_call<-list(fold=e$fold,stage=e$stage,call_index=e$gam_count,n=length(Y),family=family$family,response_range=range(Y),
    zeros=sum(Y==0),ones=sum(Y==1),weights_range=range(obsWeights))
 }),where=asNamespace('SuperLearner'),print=FALSE))
 suppressMessages(trace('qr.resid',tracer=quote({
   if(any(!is.finite(y))||any(!is.finite(qr$qr))){
    e<-get('.failure_trace',.GlobalEnv)
    e$events[[length(e$events)+1L]]<-list(function_name='qr.resid',gam=e$gam_call,
      nonfinite_y=sum(!is.finite(y)),nonfinite_qr=sum(!is.finite(qr$qr)),
      stack=vapply(sys.calls(),function(x)paste(deparse(x),collapse=' '),''))
   }
 }),where=asNamespace('stats'),print=FALSE))
 set_stream(seed)
 fit<-withCallingHandlers(capture_fit(kennedy_fit(dat,folds)),error=function(err){
  e$errors[[length(e$errors)+1L]]<-list(message=conditionMessage(err),call=deparse(conditionCall(err)),gam=e$gam_call,stack=vapply(sys.calls(),function(x)paste(deparse(x),collapse=' '),''))
 })
 suppressMessages(untrace('qr.resid',where=asNamespace('stats')));suppressMessages(untrace('SL.gam',where=asNamespace('SuperLearner')));suppressMessages(untrace('SuperLearner',where=asNamespace('SuperLearner')));suppressMessages(untrace('kennedy_fold'))
 stopifnot(is.null(fit$value),identical(fit$error,z$error))
 info<-as.list(e);atomic_save(list(key=z$key,seed=seed,reproduced_error=fit$error,trace=info),file.path('results/diagnostics',paste0(z$key,'_failure_replay.rds')))
 last<-tail(e$stages,1)[[1]];gc<-e$gam_call
 replays[[j]]<-data.frame(key=z$key,scenario=z$scenario,n=z$n,level=z$level,outer_fold=last$fold,stage=last$stage,
   training_n=last$n,observed0=last$observed0,observed1=last$observed1,last_gam_n=gc$n,
   last_gam_zero=gc$zeros,last_gam_one=gc$ones,last_gam_call_index=gc$call_index,
   nonfinite_qr_events=length(e$events),error=fit$error,
   warnings=paste(unique(fit$warnings),collapse=' | '))
}
write.csv(if(length(replays))do.call(rbind,replays)else data.frame(),'results/diagnostics/failure_replays.csv',row.names=FALSE)
writeLines('All worst-fit scores reproduced algebraically; failed Kennedy fits replayed with unchanged data and seeds.','validation/DIAGNOSTICS_VERIFIED')
