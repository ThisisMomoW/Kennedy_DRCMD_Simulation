# Aggregate saved attempts only. This script never fits an estimator.
source('R/load.R');m<-check_study()
full<-!('partial'%in%commandArgs(TRUE));rows<-diagnostics<-warnings<-list();ii<-jj<-ww<-0L
scalar<-function(z,n,default=NA_real_)if(is.null(z[[n]]))default else as.numeric(z[[n]])
text_field<-function(x)paste(unique(x),collapse=' | ')
flatten_warnings<-function(x){
 if(is.list(x)&&!is.null(x$message))return(list(x))
 if(is.list(x))return(unlist(lapply(x,flatten_warnings),recursive=FALSE))
 list()
}
for(i in seq_len(nrow(m$jobs))){
 job<-m$jobs[i,];path<-function(part)checkpoint_path(job,part)
 if(!file.exists(path('complete'))){if(full)stop('Incomplete attempt: ',job$key);next}
 s<-m$settings[job$setting,];d<-readRDS(path('data'));dat<-d$observed;f<-readRDS(path('folds'));meta<-readRDS(path('complete'))
 stopifnot(identical(object_md5(dat),meta$data_hash),identical(object_md5(f),meta$fold_hash))
 native<-readRDS(path('drcmd'));ratio<-readRDS(path('original'));direct<-readRDS(path('direct_m_gam'))
 values<-list(original=ratio$value$result,direct_m_gam=direct$result,drcmd=native)
 caps<-list(original=ratio,direct_m_gam=direct,drcmd=native$fit)
 for(mm in method_ids){
  v<-values[[mm]];cap<-caps[[mm]];estimate<-scalar(v,'estimate');se<-scalar(v,'se');lower<-scalar(v,'lower');upper<-scalar(v,'upper')
  point<-is.finite(estimate);ci<-point&&all(is.finite(c(se,lower,upper)))&&se>=0&&lower<=upper
  reason<-if(point)if(ci)''else'interval_unavailable' else if(mm=='direct_m_gam')direct$failure_type else if(nzchar(cap$error))'fitting_error'else'nonfinite_estimate'
  seed<-meta$seeds[[if(mm=='drcmd')'drcmd'else'kennedy']]
  ii<-ii+1L;rows[[ii]]<-data.frame(key=job$key,scenario=s$scenario,group=s$group,n=s$n,level=s$level,delta=s$delta,q=s$q,truth=s$truth,rep=job$rep,
   method=mm,estimate=estimate,se=se,lower=lower,upper=upper,point_valid=point,interval_valid=ci,
   realized_missingness=mean(dat$R==0),observed0=sum(dat$R==1&dat$A_obs==0,na.rm=TRUE),observed1=sum(dat$R==1&dat$A_obs==1,na.rm=TRUE),
   failure_type=reason,error=text_field(cap$error),warnings=text_field(cap$warnings),seconds=scalar(cap,'seconds'),
   data_hash=meta$data_hash,fold_hash=meta$fold_hash,generation_seed=paste(meta$seeds$data,collapse=' '),fitting_seed=paste(seed,collapse=' '))
  details<-if(mm=='direct_m_gam')unlist(lapply(direct$fold_details,function(z)lapply(z$fits,`[[`,'warning_details')),recursive=FALSE)else cap$warning_details
  details<-flatten_warnings(details)
  if(length(details))for(w in details){
   ww<-ww+1L;warnings[[ww]]<-data.frame(key=job$key,method=mm,fold=if(!length(w$fold))NA_integer_ else w$fold,
     stage=if(!length(w$stage))''else w$stage,role=if(!length(w$role))''else w$role,message=w$message)
  }
 }
 for(k in seq_along(f)){
  tr<-dat[f[[k]]$train,];base<-data.frame(key=job$key,scenario=s$scenario,group=s$group,n=s$n,level=s$level,rep=job$rep,fold=k,
    training_n=nrow(tr),evaluation_n=length(f[[k]]$test),observed0=sum(tr$R==1&tr$A_obs==0,na.rm=TRUE),observed1=sum(tr$R==1&tr$A_obs==1,na.rm=TRUE))
  if(!is.null(ratio$value)){
   p<-ratio$value$fold_details[[k]]$predictions;ds<-direct$fold_details[[k]]
   jj<-jj+1L;diagnostics[[jj]]<-cbind(base,data.frame(kappa_min=min(p$kappa_raw),kappa_max=max(p$kappa_raw),kappa_truncated=sum(p$kappa_raw!=p$kappa),
    gamma_min=min(p$gamma1_raw),gamma_max=max(p$gamma1_raw),gamma_truncated=sum(p$gamma1_raw!=p$gamma1),
    ratio_mean_min=min(p$m0,p$m1),ratio_mean_max=max(p$m0,p$m1),lambda_min=min(p$lambda1),lambda_max=max(p$lambda1),
    score_abs_max=max(abs(ratio$value$scores[f[[k]]$test,])),
    direct_converged=if(is.null(ds))NA else all(vapply(ds$fits,function(x)isTRUE(x$value$convergence$criterion_met),TRUE))))
  }else{
   jj<-jj+1L;diagnostics[[jj]]<-cbind(base,data.frame(kappa_min=NA,kappa_max=NA,kappa_truncated=NA,gamma_min=NA,gamma_max=NA,gamma_truncated=NA,
     ratio_mean_min=NA,ratio_mean_max=NA,lambda_min=NA,lambda_max=NA,score_abs_max=NA,direct_converged=NA))
  }
 }
}
r<-do.call(rbind,rows);stopifnot(!anyDuplicated(r[c('key','method')]))
if(full)stopifnot(nrow(r)==9600,length(unique(r$key))==3200)
summaries<-lapply(split(r,interaction(r$scenario,r$method,drop=TRUE)),function(z){
 p<-z[z$point_valid,];ci<-z[z$interval_valid,];err<-p$estimate-p$truth;B<-length(err);Bci<-nrow(ci)
 rmse<-sqrt(mean(err^2));coverage<-mean(ci$lower<=ci$truth&ci$truth<=ci$upper);ss<-sort(err^2,decreasing=TRUE)
 data.frame(z[1,c('scenario','group','n','level','delta','q','truth','method')],attempted=nrow(z),point_n=B,ci_n=Bci,
  point_failures=nrow(z)-B,interval_failures=nrow(z)-Bci,rd_outside=sum(abs(p$estimate)>1),
  bias=mean(err),bias_mcse=sd(err)/sqrt(B),rmse=rmse,rmse_mcse=sd(err^2)/(2*rmse*sqrt(B)),
  coverage=coverage,coverage_mcse=sqrt(coverage*(1-coverage)/Bci),
  minimum=min(p$estimate),maximum=max(p$estimate),top1_squared_error_fraction=ss[1]/sum(ss),top5_squared_error_fraction=sum(head(ss,5))/sum(ss),
  mean_missingness=mean(z$realized_missingness),mean_seconds=mean(z$seconds))
})
s<-do.call(rbind,summaries);rownames(s)<-NULL;s<-s[order(s$group,s$level,s$n,match(s$method,method_ids)),]
if(full)stopifnot(nrow(s)==48,all(s$attempted==200))
paired<-list();q<-0L
for(id in unique(r$scenario))for(mm in c('original','direct_m_gam')){
 a<-subset(r,scenario==id&method=='drcmd'&point_valid);b<-subset(r,scenario==id&method==mm&point_valid)
 z<-merge(a,b,by='key');delta<-(z$estimate.x-z$truth.x)^2-(z$estimate.y-z$truth.y)^2
 q<-q+1L;paired[[q]]<-data.frame(scenario=id,group=a$group[1],n=a$n[1],level=a$level[1],comparator=mm,paired_n=nrow(z),
   difference=mean(delta),mcse=sd(delta)/sqrt(length(delta)),comparator_smaller_abs_error=mean(abs(z$estimate.y-z$truth.y)<abs(z$estimate.x-z$truth.x)))
}
paired<-do.call(rbind,paired);diag<-do.call(rbind,diagnostics)
warn<-if(length(warnings))do.call(rbind,warnings)else data.frame(key=character(),method=character(),fold=integer(),stage=character(),role=character(),message=character())
ref<-read.csv('reference/summary.csv');cmp<-merge(s,ref,by=c('scenario','method'),suffixes=c('_new','_reference'))
for(metric in c('bias','rmse','coverage')){
 cmp[[paste0(metric,'_change')]]<-cmp[[paste0(metric,'_new')]]-cmp[[paste0(metric,'_reference')]]
 cmp[[paste0(metric,'_change_mcse')]]<-sqrt(cmp[[paste0(metric,'_mcse_new')]]^2+cmp[[paste0(metric,'_mcse_reference')]]^2)
}
write.csv(r,'results/replicates.csv',row.names=FALSE);write.csv(s,'results/summary.csv',row.names=FALSE)
write.csv(paired,'results/paired_mse.csv',row.names=FALSE);write.csv(cmp,'results/reference_comparison.csv',row.names=FALSE)
write.csv(diag,'results/fold_diagnostics.csv',row.names=FALSE);write.csv(warn,'results/warnings.csv',row.names=FALSE)
write.csv(subset(r,!point_valid|!interval_valid),'results/failures.csv',row.names=FALSE)
extremes<-do.call(rbind,lapply(split(subset(r,point_valid),interaction(subset(r,point_valid)$scenario,subset(r,point_valid)$method,drop=TRUE)),function(z){
 z$squared_error<-(z$estimate-z$truth)^2;z[head(order(z$squared_error,decreasing=TRUE),5),]
}))
write.csv(extremes,'results/largest_errors.csv',row.names=FALSE)
for(g in 1:3){path<-file.path(scenario_paths[g],'results')
 write.csv(subset(r,group==g),file.path(path,'replicates.csv'),row.names=FALSE)
 write.csv(subset(s,group==g),file.path(path,'summary.csv'),row.names=FALSE)
 write.csv(subset(paired,group==g),file.path(path,'paired_mse.csv'),row.names=FALSE)
 write.csv(subset(diag,group==g),file.path(path,'fold_diagnostics.csv'),row.names=FALSE)
}
if(full)writeLines(c('3200 distinct attempted datasets; 9600 method records; 48 summary rows.',
 'All finite estimates retained. Coverage uses valid intervals. Data/fold hashes agree.'),'validation/AGGREGATION_VERIFIED')
cat('Aggregated',length(unique(r$key)),'datasets and',nrow(s),'summary rows.\n')
