eem_package_call<-function(fun){
 found<-list()
 visit<-function(x){if(missing(x))return(invisible(NULL));if(is.call(x)){
  if(identical(x[[1]],quote(drcmd::drcmd)))found[[length(found)+1L]]<<-x
  for(y in as.list(x))visit(y)
 }}
 visit(body(fun));stopifnot(length(found)==1);found[[1]]
}
eem_verify<-function(){
 x<-eem_initialize();b<-check_study()
 # Enforce the requested one-argument change directly on the wrapper call AST.
 ordinary<-eem_package_call(drcmd_point);expected<-ordinary;expected$eem_ind<-TRUE
 stopifnot(identical(expected,eem_package_call(drcmd_eem_point)),identical(ordinary$eem_ind,FALSE))
 ids<-which(!duplicated(b$jobs$group));stopifnot(length(ids)==3)
 rows<-list();values<-list();baseline_csv<-read.csv('results/replicates.csv')
 eem_start_observers()
 for(j in seq_along(ids)){
  i<-ids[j];input<-eem_read_input(i,b);job<-input$job
  checkpoint_files<-vapply(c('original','direct_m_gam','drcmd','data','folds','complete'),function(p)checkpoint_path(job,p),'')
  before<-tools::md5sum(checkpoint_files)
  out<-eem_one(i,b,x);after<-tools::md5sum(checkpoint_files)
  stopifnot(identical(before,after),isTRUE(out$fold_match),out$record$point_valid,out$record$interval_valid,
    out$fit$audit$counts['est_varphi_eem']==5,out$fit$audit$counts['est_varphi']==0,out$fit$audit$counts['est_psi_tml']==0)
  cap<-out$fit$fit$value;params<-cap$params
  stopifnot(identical(params$eem_ind,TRUE),identical(params$tml,FALSE),params$k==5,
    all(is.na(params$Rprobs)),params$default_learners=='SL.gam')
  fits<-out$fit$audit$sl;aug<-Filter(function(z)z$stage=='est_varphi_eem',fits)
  stopifnot(length(aug)==15,all(vapply(aug,function(z)z$family=='gaussian'&&identical(z$library,'SL.gam')&&z$requested_V==5,TRUE)))
  stopifnot(identical(vapply(aug,`[[`,'','role'),rep(c('arm1','arm0','unused_direct_contrast'),5)))
  # The two arms use the EEM weights (R/kappa-1)^2 and all training rows.
  for(k in 1:5){
   raw<-Filter(function(z)z$fold==k&&z$target=='truncate_r',out$fit$audit$raw)
   stopifnot(length(raw)==1)
   kappa<-pmin(.99,pmax(.01,raw[[1]]$raw));idx<-input$folds[[k]]$train
   weight<-(input$dat$R/kappa-1)^2
   these<-Filter(function(z)z$fold==k,aug)
   for(v in these)close_numeric(v$weights,weight[idx])
  }
  close_numeric(out$record$estimate,as.numeric(cap$results$estimates$psi_hat_ate))
  close_numeric(out$record$se,as.numeric(cap$results$ses$psi_hat_ate))
  close_numeric(out$record$lower,out$record$estimate-1.96*out$record$se)
  close_numeric(out$record$upper,out$record$estimate+1.96*out$record$se)
  # Read the three untouched checkpoint estimates; do not refit them for this extension.
  original<-readRDS(checkpoint_path(job,'original'))$value$result
  direct<-readRDS(checkpoint_path(job,'direct_m_gam'))$result
  dr<-readRDS(checkpoint_path(job,'drcmd'))
  vals<-list(original=original,direct_m_gam=direct,drcmd=dr)
  max_error<-0
  for(mm in names(vals)){
   ref<-subset(baseline_csv,key==job$key&method==mm)
   actual<-unlist(vals[[mm]][c('estimate','se','lower','upper')]);wanted<-as.numeric(ref[,c('estimate','se','lower','upper')])
   close_numeric(as.numeric(actual),wanted);max_error<-max(max_error,max(abs(actual-wanted)))
   values[[length(values)+1L]]<-data.frame(key=job$key,method=mm,estimate=actual[1],se=actual[2],lower=actual[3],upper=actual[4])
  }
  values[[length(values)+1L]]<-out$record[,c('key','method','estimate','se','lower','upper')]
  rows[[j]]<-data.frame(key=job$key,group=job$group,data_match=TRUE,seed_match=TRUE,fold_match=TRUE,
    native_eem_calls=5,ordinary_augmentation_calls=0,tml_calls=0,original_checkpoint_md5_unchanged=TRUE,
    original_estimate_se_ci_max_error=max_error)
 }
 # Meaningful aggregation check: unavailable intervals and nonfinite points have separate denominators.
 z<-out$record[rep(1,3),];z$estimate<-c(z$truth[1]+1,z$truth[1]+3,NA_real_)
 z$point_valid<-c(TRUE,TRUE,FALSE);z$interval_valid<-c(TRUE,FALSE,FALSE)
 z$lower<-z$truth-1;z$upper<-z$truth+1
 ss<-eem_summary_row(z,b$settings[1,],'drcmd_eem',3L)
 close_numeric(ss$bias,2);close_numeric(ss$rmse,sqrt(5));stopifnot(ss$ci_n==1,ss$coverage==1,ss$point_failures==1,ss$interval_failures==2)
 compare<-eem_summarize(partial=TRUE)
 stopifnot(setequal(unique(compare$replicates$method),c(method_ids,'drcmd_eem')),nrow(compare$summary)==64)
 for(i in ids){selected_key<-b$jobs$key[i];z<-compare$replicates[compare$replicates$key==selected_key,];stopifnot(nrow(z)==4)
  a<-subset(z,method=='drcmd');d<-subset(z,method=='drcmd_eem')
  ps<-subset(compare$paired,scenario==a$scenario)
  if(ps$paired_n==1)close_numeric(ps$difference,(a$estimate-a$truth)^2-(d$estimate-d$truth)^2)
 }
 write.csv(do.call(rbind,rows),eem_path('verification/checks.csv'),row.names=FALSE)
 write.csv(do.call(rbind,values),eem_path('verification/four_method_examples.csv'),row.names=FALSE)
 writeLines(c('PASS: only eem_ind changes in the exported package call.',
  'PASS: one fixed existing replicate per scenario, original data/seed/folds match.',
  'PASS: 5 EEM calls and 0 ordinary augmentation/TMLE calls per verified fit.',
  'PASS: Gaussian GAM EEM fits use the native weights and internal CV settings.',
  'PASS: original three checkpoint hashes and all estimate/SE/CI values preserved.',
  'PASS: four methods in separate outputs; failure denominators and paired MSE checked.',
  'Only the three verification EEM attempts were requested here; no full jobs submitted.'),eem_path('verification/VERIFIED'))
 eem_guard();cat(readLines(eem_path('verification/VERIFIED')),sep='\n')
}
