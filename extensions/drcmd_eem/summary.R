eem_mean<-function(x)if(length(x))mean(x)else NA_real_
eem_summary_row<-function(z,setting,method,planned=200L) {
 p<-z[z$point_valid,];ci<-z[z$interval_valid,];err<-p$estimate-p$truth;B<-length(err);Nci<-nrow(ci)
 rmse<-sqrt(eem_mean(err^2));coverage<-eem_mean(ci$lower<=ci$truth&ci$truth<=ci$upper)
 ss<-sort(err^2,decreasing=TRUE)
 data.frame(setting[,c('scenario','group','n','level','delta','q','truth')],method=method,
  attempted=nrow(z),point_n=B,ci_n=Nci,point_failures=nrow(z)-B,interval_failures=nrow(z)-Nci,rd_outside=sum(abs(p$estimate)>1),
  bias=eem_mean(err),bias_mcse=if(B>1)sd(err)/sqrt(B)else NA_real_,rmse=rmse,
  rmse_mcse=if(B>1&&rmse>0)sd(err^2)/(2*rmse*sqrt(B))else NA_real_,
  coverage=coverage,coverage_mcse=if(Nci)sqrt(coverage*(1-coverage)/Nci)else NA_real_,
  minimum=if(B)min(p$estimate)else NA_real_,maximum=if(B)max(p$estimate)else NA_real_,
  top1_squared_error_fraction=if(B&&sum(ss)>0)ss[1]/sum(ss)else NA_real_,
  top5_squared_error_fraction=if(B&&sum(ss)>0)sum(head(ss,5))/sum(ss)else NA_real_,
  mean_missingness=eem_mean(z$realized_missingness),mean_seconds=eem_mean(z$seconds),
  planned=planned,pending=planned-nrow(z),fit_failures=sum(!z$point_valid),
  point_failure_rate=if(nrow(z))sum(!z$point_valid)/nrow(z)else NA_real_,
  interval_failure_rate=if(nrow(z))sum(!z$interval_valid)/nrow(z)else NA_real_)
}
eem_paired<-function(r,settings) {
 records<-summaries<-list()
 for(i in seq_len(nrow(settings))){
  s<-settings[i,];a<-subset(r,scenario==s$scenario&method=='drcmd'&point_valid)
  b<-subset(r,scenario==s$scenario&method=='drcmd_eem'&point_valid)
  z<-merge(a,b,by='key',suffixes=c('_drcmd','_eem'))
  stopifnot(identical(z$data_hash_drcmd,z$data_hash_eem),identical(z$fold_hash_drcmd,z$fold_hash_eem),
    identical(z$fitting_seed_drcmd,z$fitting_seed_eem))
  d<-(z$estimate_drcmd-z$truth_drcmd)^2-(z$estimate_eem-z$truth_eem)^2
  records[[i]]<-data.frame(key=z$key,scenario=rep(s$scenario,nrow(z)),estimate_drcmd=z$estimate_drcmd,
    estimate_eem=z$estimate_eem,truth=z$truth_drcmd,squared_error_difference=d)
  summaries[[i]]<-data.frame(s[,c('scenario','group','n','level','delta')],paired_n=nrow(z),
    difference=eem_mean(d),mcse=if(length(d)>1)sd(d)/sqrt(length(d))else NA_real_,
    eem_smaller_absolute_error=eem_mean(abs(z$estimate_eem-z$truth_eem)<abs(z$estimate_drcmd-z$truth_drcmd)))
 }
 list(summary=do.call(rbind,summaries),replicates=do.call(rbind,records))
}
eem_summarize<-function(partial=FALSE) {
 x<-eem_guard();b<-check_study()
 original<-read.csv('results/replicates.csv');baseline<-read.csv('results/summary.csv')
 rows<-warn<-checks<-list()
 for(i in seq_len(nrow(x$jobs))){
  job<-x$jobs[i,];path<-eem_path('checkpoints',paste0(job$key,'.rds'))
  if(!file.exists(path))next
  v<-readRDS(path);r<-v$record;ref<-subset(original,key==job$key&method=='drcmd')
  stopifnot(nrow(r)==1,nrow(ref)==1,identical(v$source_pin,x$sources),isTRUE(v$fold_match),
    identical(r$data_hash,ref$data_hash),identical(r$fold_hash,ref$fold_hash),identical(r$fitting_seed,ref$fitting_seed),
    identical(r$generation_seed,ref$generation_seed),identical(v$seed,b$seeds[[i]]$drcmd),
    v$fit$audit$counts['drcmd']==1,v$fit$audit$counts['est_varphi']==0,v$fit$audit$counts['est_psi_tml']==0)
  stopifnot(identical(v$input_hashes,tools::md5sum(names(v$input_hashes))))
  if(!is.null(v$fit$fit$value))stopifnot(v$fit$audit$counts['est_varphi_eem']==5,isTRUE(v$fit$fit$value$params$eem_ind))
  rows[[length(rows)+1L]]<-r
  for(w in v$fit$fit$warning_details)warn[[length(warn)+1L]]<-data.frame(key=job$key,method='drcmd_eem',fold=w$fold,stage=w$stage,role=w$role,message=w$message)
  checks[[length(checks)+1L]]<-data.frame(key=job$key,data_match=TRUE,fold_match=TRUE,seed_match=TRUE,
    eem_calls=unname(v$fit$audit$counts['est_varphi_eem']),standard_augmentation_calls=unname(v$fit$audit$counts['est_varphi']))
 }
 eem<-if(length(rows))do.call(rbind,rows)else original[FALSE,]
 if(!partial&&nrow(eem)!=3200)stop('EEM has ',nrow(eem),' of 3200 attempts. Use summarize --partial for explicitly incomplete tables.')
 stopifnot(!anyDuplicated(eem$key),identical(names(eem),names(original)))
 r<-rbind(original,eem);stopifnot(!anyDuplicated(r[c('key','method')]))
 methods<-c(method_ids,'drcmd_eem');summary<-list();matched_summary<-list()
 matched<-r[r$key%in%eem$key,]
 for(i in seq_len(nrow(x$settings)))for(mm in methods){
  s<-x$settings[i,];z<-r[r$scenario==s$scenario&r$method==mm,]
  summary[[length(summary)+1L]]<-eem_summary_row(z,s,mm)
  mz<-matched[matched$scenario==s$scenario&matched$method==mm,]
  if(nrow(mz))matched_summary[[length(matched_summary)+1L]]<-eem_summary_row(mz,s,mm,planned=sum(eem$scenario==s$scenario))
 }
 summary<-do.call(rbind,summary)
 # Recompute the three original methods and require agreement with every saved summary field.
 actual<-summary[summary$method!='drcmd_eem',]
 expected<-baseline[match(paste(actual$scenario,actual$method),paste(baseline$scenario,baseline$method)),]
 for(n in names(baseline)){
  if(is.numeric(expected[[n]]))close_numeric(actual[[n]],expected[[n]])else stopifnot(identical(actual[[n]],expected[[n]]))
 }
 paired<-eem_paired(r,x$settings);out<-eem_path('comparison');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 write.csv(eem,eem_path('eem_replicates.csv'),row.names=FALSE)
 write.csv(r,file.path(out,'replicates.csv'),row.names=FALSE)
 write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
 write.csv(matched,file.path(out,'matched_replicates.csv'),row.names=FALSE)
 write.csv(if(length(matched_summary))do.call(rbind,matched_summary)else summary[FALSE,],file.path(out,'matched_summary.csv'),row.names=FALSE)
 write.csv(paired$summary,file.path(out,'paired_mse.csv'),row.names=FALSE)
 write.csv(paired$replicates,file.path(out,'paired_squared_errors.csv'),row.names=FALSE)
 write.csv(subset(eem,!point_valid|!interval_valid),file.path(out,'failures.csv'),row.names=FALSE)
 write.csv(if(length(warn))do.call(rbind,warn)else data.frame(key=character(),method=character(),fold=integer(),stage=character(),role=character(),message=character()),file.path(out,'warnings.csv'),row.names=FALSE)
 write.csv(if(length(checks))do.call(rbind,checks)else data.frame(),file.path(out,'identity_checks.csv'),row.names=FALSE)
 for(g in 1:3){p<-file.path(out,basename(scenario_paths[g]));dir.create(p,showWarnings=FALSE)
  write.csv(subset(r,group==g),file.path(p,'replicates.csv'),row.names=FALSE)
  write.csv(subset(summary,group==g),file.path(p,'summary.csv'),row.names=FALSE)
  write.csv(subset(paired$summary,group==g),file.path(p,'paired_mse.csv'),row.names=FALSE)
 }
 writeLines(c(if(nrow(eem)==3200)'COMPLETE'else'PARTIAL',paste('EEM attempts:',nrow(eem),'of 3200'),
  'Original 9600 method records and 48 summaries preserved.',
  'Paired difference = squared error(drcmd) - squared error(drcmd_eem); positive favors EEM.',
  'Pending EEM fits are not counted as failures. Coverage uses valid intervals.'),file.path(out,'STATUS.txt'))
 eem_guard();cat('Saved four-method tables:',nrow(eem),'EEM attempts,',nrow(summary),'planned setting-method rows.\n')
 invisible(list(replicates=r,summary=summary,paired=paired$summary))
}
