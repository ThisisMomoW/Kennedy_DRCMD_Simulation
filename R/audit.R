# Read-only tracing. Estimator bodies and installed packages are unchanged.
# Context capture adds warning provenance without changing warning handling or RNG.
study_context <- new.env(parent=emptyenv());study_context$method<-'';study_context$fold<-0L
capture_fit <- function(expr) {
 w<-msg<-character();details<-list();err<-'';t<-proc.time()[3];value<-NULL
 console<-capture.output({value<-tryCatch(withCallingHandlers(expr,
  warning=function(x) {
   message<-conditionMessage(x);w<<-c(w,message)
   e<-if(exists('.gam_audit',.GlobalEnv))get('.gam_audit',.GlobalEnv) else NULL
   native<-study_context$method=='drcmd'&&!is.null(e)
   details[[length(details)+1L]]<<-list(method=study_context$method,
    fold=if(native)e$fold else study_context$fold,stage=if(native)e$stage else 'native_or_local_fit',
    role=if(native)e$role else '',message=message)
   invokeRestart('muffleWarning')
  },message=function(x){msg<<-c(msg,conditionMessage(x));invokeRestart('muffleMessage')}),
  error=function(x){err<<-conditionMessage(x);NULL})})
 list(value=value,warnings=w,warning_details=details,messages=msg,error=err,console=console,seconds=unname(proc.time()[3]-t))
}
study_audit_reset <- function() {
 assign('.gam_audit',new.env(parent=emptyenv()),.GlobalEnv);reset_audit()
 e<-get('.gam_audit',.GlobalEnv);e$active<-FALSE;e$fold<-0L;e$role<-'';e$augmentation<-0L
 e$raw<-e$centering<-e$fold_results<-list()
}
study_start_audit <- function() {
 study_audit_reset()
 for(n in c('drcmd','est_kappa','est_g','est_m_a','est_varphi','est_varphi_eem','est_psi_tml')) {
  expr<-substitute({e<-get('.gam_audit',.GlobalEnv);if(e$active){e$counts[N]<-e$counts[N]+1L;e$stage<-N;e$role<-'';if(N=='est_varphi')e$augmentation<-0L}},list(N=n))
  suppressMessages(trace(n,tracer=expr,where=asNamespace('drcmd'),print=FALSE))
 }
 suppressMessages(trace('create_folds',exit=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active)e$folds<-returnValue()}),where=asNamespace('drcmd'),print=FALSE))
 suppressMessages(trace('drcmd_est_fold',tracer=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active)e$fold<-e$fold+1L}),
  exit=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active)e$fold_results[[e$fold]]<-returnValue()}),where=asNamespace('drcmd'),print=FALSE))
 suppressMessages(trace('get_phi_hat',exit=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active){v<-returnValue();e$centering[[e$fold]]<-list(plugin0=v$plugin0,plugin1=v$plugin1,training_ids=idx)}}),where=asNamespace('drcmd'),print=FALSE))
 suppressMessages(trace('fit_augmentation_sl',tracer=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active){e$augmentation<-e$augmentation+1L;e$role<-c('unused_direct_contrast','arm1','arm0')[e$augmentation]}}),where=asNamespace('drcmd'),print=FALSE))
 for(n in c('truncate_r','truncate_g')) {
  expr<-substitute({e<-get('.gam_audit',.GlobalEnv);if(e$active){
   e$raw[[length(e$raw)+1L]]<-list(fold=e$fold,target=N,raw=as.numeric(x),cutoff=cutoff)
   e$clipping[[length(e$clipping)+1L]]<-list(fold=e$fold,function_name=N,n=length(x),clipped=sum(x<cutoff|x>1-cutoff),min=min(x),max=max(x),cutoff=cutoff)
  }},list(N=n))
  suppressMessages(trace(n,tracer=expr,where=asNamespace('drcmd'),print=FALSE))
 }
 suppressMessages(trace('SuperLearner',tracer=quote({
  e<-get('.gam_audit',.GlobalEnv);if(e$active){ww<-obsWeights;if(is.null(ww))ww<-rep(1,length(Y))
   e$sl[[length(e$sl)+1L]]<-list(fold=e$fold,stage=e$stage,role=e$role,n=length(Y),
    weights=as.numeric(ww),response=as.numeric(Y),family=family$family,library=SL.library,requested_V=cvControl$V)
  }
 }),where=asNamespace('SuperLearner'),print=FALSE))
 suppressMessages(trace('CVFolds',exit=quote({e<-get('.gam_audit',.GlobalEnv);if(e$active){f<-returnValue();
  e$cv[[length(e$cv)+1L]]<-list(fold=e$fold,stage=e$stage,role=e$role,n=N,requested_V=cvControl$V,actual_V=length(f),evaluation_counts=lengths(f),training_counts=N-lengths(f))
 }}),where=asNamespace('SuperLearner'),print=FALSE))
 suppressMessages(trace('kennedy_fold',tracer=quote({study_context$fold<-study_context$fold+1L}),where=.GlobalEnv,print=FALSE))
}
