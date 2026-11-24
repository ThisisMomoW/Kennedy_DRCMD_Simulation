# Separate outputs and source pins; original R/ and its manifest remain immutable.
eem_root<-function()'results/drcmd_eem_v1'
eem_path<-function(...)file.path(eem_root(),...)
eem_sources<-function() {
 f<-c('run.R',list.files('extensions/drcmd_eem',pattern='\\.R$',full.names=TRUE),
      'scripts/eem/slurm.sh','scripts/eem/submit.sh')
 tools::md5sum(sort(f))
}
eem_protected<-function() {
 c(names(source_fingerprints()),'results/manifest.rds','results/attempt_plan.csv',
   list.files('results',pattern='\\.csv$',full.names=TRUE),
   list.files('reports',pattern='\\.(Rmd|tex|pdf)$',full.names=TRUE),'reports/support.R',
   file.path(scenario_paths,'settings.csv'),
   unlist(lapply(file.path(scenario_paths,'results'),list.files,pattern='\\.csv$',full.names=TRUE)))
}
eem_initialize<-function() {
 if(file.exists(eem_path('manifest.rds')))return(eem_guard())
 b<-check_study();r<-read.csv('results/replicates.csv');s<-read.csv('results/summary.csv')
 stopifnot(nrow(r)==9600,nrow(s)==48,all(s$attempted==200),setequal(unique(r$method),method_ids))
 for(p in c('checkpoints','logs','verification','comparison','report'))dir.create(eem_path(p),recursive=TRUE,showWarnings=FALSE)
 protected<-unique(eem_protected());stopifnot(all(file.exists(protected)))
 x<-list(id='drcmd_eem_v1',method='drcmd_eem',base_manifest_md5=unname(tools::md5sum('results/manifest.rds')),
   sources=eem_sources(),protected=tools::md5sum(protected),software=software_info(),R=R.version.string,
   jobs=b$jobs,settings=b$settings,created=as.character(Sys.time()),
   arguments=list(default_learners='SL.gam',eem_ind=TRUE,tml=FALSE,k=5,cv_folds=5,cutoff=.01,parallel=FALSE),
   estimate='results$estimates$psi_hat_ate',se='results$ses$psi_hat_ate',interval='estimate +/- 1.96*SE')
 atomic_save(x,eem_path('manifest.rds'))
 write.csv(x$software,eem_path('software.csv'),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),eem_path('session_info.txt'))
 x
}
eem_guard<-function() {
 b<-check_study();x<-readRDS(eem_path('manifest.rds'))
 stopifnot(identical(x$sources,eem_sources()),
   identical(x$protected,tools::md5sum(names(x$protected))),
   identical(x$base_manifest_md5,unname(tools::md5sum('results/manifest.rds'))),
   identical(x$software[,1:3],software_info()[,1:3]),identical(x$R,R.version.string),
   identical(x$jobs,b$jobs),identical(x$settings,b$settings))
 x
}
eem_read_input<-function(index,b=check_study()) {
 job<-b$jobs[index,];path<-function(part)checkpoint_path(job,part)
 files<-setNames(vapply(c('data','folds','drcmd','complete'),path,''),c('data','folds','drcmd','complete'))
 stopifnot(all(file.exists(files)));hashes<-tools::md5sum(unname(files))
 d<-readRDS(files['data']);folds<-readRDS(files['folds']);dr<-readRDS(files['drcmd']);meta<-readRDS(files['complete'])
 verify_folds(folds,nrow(d$observed))
 stopifnot(identical(d$seed,b$seeds[[index]]$data),identical(meta$seeds,b$seeds[[index]]),
   identical(meta$data_hash,d$hash),identical(d$hash,object_md5(d$observed)),
   identical(meta$fold_hash,object_md5(folds)),identical(folds,dr$audit$folds),
   identical(is.na(d$observed$A_obs),d$observed$R==0))
 list(job=job,setting=b$settings[job$setting,],dat=d$observed,seed=meta$seeds$drcmd,
   generation_seed=d$seed,folds=folds,data_hash=d$hash,fold_hash=meta$fold_hash,
   baseline_drcmd=dr,files=files,hashes=hashes)
}
eem_start_observers<-function() {
 study_start_audit()
 # Native non-EEM order: contrast, arm1, arm0. Native EEM order: arm1, arm0, contrast.
 suppressMessages(trace('est_varphi_eem',tracer=quote({
  e<-get('.gam_audit',.GlobalEnv)
  if(e$active){e$counts['est_varphi_eem']<-e$counts['est_varphi_eem']+1L;e$stage<-'est_varphi_eem';e$role<-'';e$augmentation<-0L}
 }),where=asNamespace('drcmd'),print=FALSE))
 suppressMessages(trace('fit_augmentation_sl',tracer=quote({
  e<-get('.gam_audit',.GlobalEnv)
  if(e$active){e$augmentation<-e$augmentation+1L
   roles<-if(e$stage=='est_varphi_eem')c('arm1','arm0','unused_direct_contrast')else c('unused_direct_contrast','arm1','arm0')
   e$role<-roles[e$augmentation]
  }
 }),where=asNamespace('drcmd'),print=FALSE))
}
eem_numeric<-function(z,n)if(is.null(z[[n]]))NA_real_ else as.numeric(z[[n]])
eem_record<-function(input,fit) {
 z<-input;d<-z$dat;s<-z$setting
 est<-eem_numeric(fit,'estimate');se<-eem_numeric(fit,'se');lo<-eem_numeric(fit,'lower');hi<-eem_numeric(fit,'upper')
 point<-is.finite(est);ci<-point&&all(is.finite(c(se,lo,hi)))&&se>=0&&lo<=hi
 error<-paste(unique(fit$fit$error),collapse=' | ')
 failure<-if(point)if(ci)''else'interval_unavailable'else if(nzchar(error))'fitting_error'else'nonfinite_estimate'
 data.frame(key=z$job$key,scenario=s$scenario,group=s$group,n=s$n,level=s$level,delta=s$delta,q=s$q,truth=s$truth,rep=z$job$rep,
   method='drcmd_eem',estimate=est,se=se,lower=lo,upper=hi,point_valid=point,interval_valid=ci,
   realized_missingness=mean(d$R==0),observed0=sum(d$R==1&d$A_obs==0,na.rm=TRUE),observed1=sum(d$R==1&d$A_obs==1,na.rm=TRUE),
   failure_type=failure,error=error,warnings=paste(unique(fit$fit$warnings),collapse=' | '),seconds=fit$fit$seconds,
   data_hash=z$data_hash,fold_hash=z$fold_hash,generation_seed=paste(z$generation_seed,collapse=' '),fitting_seed=paste(z$seed,collapse=' '))
}
eem_one<-function(index,b=check_study(),manifest=eem_guard()) {
 input<-eem_read_input(index,b);path<-eem_path('checkpoints',paste0(input$job$key,'.rds'))
 if(file.exists(path)){
  old<-readRDS(path)
  stopifnot(identical(old$input_hashes,input$hashes),identical(old$seed,input$seed),
    identical(old$data_hash,input$data_hash),identical(old$fold_hash,input$fold_hash),identical(old$source_pin,manifest$sources))
  return(invisible(old))
 }
 log<-eem_path('logs',paste0(input$job$key,'.log'))
 cat(format(Sys.time()),'start',Sys.getenv('SLURM_JOB_ID'),'\n',file=log,append=TRUE)
 study_context$method<-'drcmd';study_context$fold<-0L;study_audit_reset()
 e<-get('.gam_audit',.GlobalEnv);e$active<-TRUE
 on.exit({e$active<-FALSE},add=TRUE)
 # Reset to the original DR-CMD stream immediately before the exported package call.
 set_stream(input$seed);fit<-drcmd_eem_point(input$dat);e$active<-FALSE
 for(j in seq_along(fit$fit$warning_details))fit$fit$warning_details[[j]]$method<-'drcmd_eem'
 record<-eem_record(input,fit)
 fold_match<-identical(fit$audit$folds,input$folds)
 stopifnot(identical(input$hashes,tools::md5sum(unname(input$files))))
 # Always checkpoint an attempted fit before checking native branch evidence.
 out<-list(id=manifest$id,record=record,fit=fit,seed=input$seed,generation_seed=input$generation_seed,
   data_hash=input$data_hash,fold_hash=input$fold_hash,folds=input$folds,fold_match=fold_match,
   input_hashes=input$hashes,source_pin=manifest$sources,software=manifest$software,
   host=Sys.info()[['nodename']],job=Sys.getenv('SLURM_JOB_ID'),time=as.character(Sys.time()))
 atomic_save(out,path)
 stopifnot(fold_match,fit$audit$counts['drcmd']==1,fit$audit$counts['est_varphi']==0,fit$audit$counts['est_psi_tml']==0)
 if(!is.null(fit$fit$value))stopifnot(isTRUE(fit$fit$value$params$eem_ind),identical(fit$fit$value$params$tml,FALSE),fit$audit$counts['est_varphi_eem']==5)
 cat(format(Sys.time()),'finished',record$failure_type,'\n',file=log,append=TRUE)
 cat(input$job$key,'EEM',if(record$point_valid)'finite'else record$failure_type,'\n');flush.console()
 invisible(out)
}
