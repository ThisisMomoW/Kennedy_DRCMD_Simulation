# Four-method saved-result validation and plots. No fitting occurs here.
report_methods<-c('original','direct_m_gam','drcmd','drcmd_eem')
report_labels<-c(original='Kennedy ratio GAM',direct_m_gam='Kennedy direct GAM',drcmd='DRCMD',drcmd_eem='DRCMD EEM')
report_colors<-c(original='#B8411A',direct_m_gam='#00867D',drcmd='#2E5378',drcmd_eem='#8E5BB7')
report_paths<-c('scenarios/01_observation_heterogeneity','scenarios/02_outcome_slope','scenarios/03_nonlinear_outcome')
report_close<-function(a,b){
 a<-as.numeric(a);b<-as.numeric(b);stopifnot(length(a)==length(b),identical(is.na(a),is.na(b)))
 keep<-!is.na(a);stopifnot(all(abs(a[keep]-b[keep])<=1e-10*(1+abs(b[keep]))))
}
read_report_data<-function(root){
 base<-file.path(root,'results/drcmd_eem_v1/comparison')
 stopifnot(readLines(file.path(base,'STATUS.txt'))[1]=='COMPLETE')
 r<-read.csv(file.path(base,'replicates.csv'));s<-read.csv(file.path(base,'summary.csv'))
 stopifnot(nrow(r)==12800,nrow(s)==64,length(unique(r$key))==3200,all(s$attempted==200),all(s$pending==0),
  setequal(unique(r$method),report_methods),!anyDuplicated(r[c('key','method')]))
 checks<-list()
 for(i in seq_len(nrow(s))){
  z<-r[r$scenario==s$scenario[i]&r$method==s$method[i],];finite<-is.finite(z$estimate)
  valid<-finite&is.finite(z$se)&z$se>=0&is.finite(z$lower)&is.finite(z$upper)&z$lower<=z$upper
  stopifnot(identical(z$point_valid,finite),identical(z$interval_valid,valid))
  e<-z$estimate[finite]-z$truth[finite];ci<-z[valid,];B<-length(e);N<-nrow(ci)
  rr<-if(B)sqrt(mean(e^2))else NA_real_;cc<-if(N)mean(ci$lower<=ci$truth&ci$truth<=ci$upper)else NA_real_
  vals<-c(bias=if(B)mean(e)else NA_real_,bias_mcse=if(B>1)sd(e)/sqrt(B)else NA_real_,
   rmse=rr,rmse_mcse=if(B>1&&rr>0)sd(e^2)/(2*rr*sqrt(B))else NA_real_,
   coverage=cc,coverage_mcse=if(N)sqrt(cc*(1-cc)/N)else NA_real_)
  report_close(vals,as.numeric(s[i,names(vals)]))
  stopifnot(nrow(z)==s$attempted[i],B==s$point_n[i],N==s$ci_n[i],sum(!finite)==s$point_failures[i],
   sum(!valid)==s$interval_failures[i],sum(abs(z$estimate[finite])>1)==s$rd_outside[i])
  checks[[i]]<-data.frame(scenario=s$scenario[i],method=s$method[i],max_difference=max(c(0,abs(vals-as.numeric(s[i,names(vals)]))),na.rm=TRUE))
 }
 for(z in split(r,r$key)){
  stopifnot(nrow(z)==4,length(unique(z$data_hash))==1,length(unique(z$fold_hash))==1)
  d<-z[z$method%in%c('drcmd','drcmd_eem'),];stopifnot(length(unique(d$fitting_seed))==1)
 }
 original<-read.csv(file.path(root,'results/replicates.csv'));old<-r[r$method!='drcmd_eem',]
 old<-old[match(paste(original$key,original$method),paste(old$key,old$method)),]
 for(n in names(original)){
  if(is.numeric(original[[n]]))report_close(old[[n]],original[[n]])else stopifnot(identical(old[[n]],original[[n]]))
 }
 paired<-read.csv(file.path(base,'paired_mse.csv'))
 for(i in seq_len(nrow(paired))){
  z<-r[r$scenario==paired$scenario[i]&r$point_valid,];a<-z[z$method=='drcmd',];b<-z[z$method=='drcmd_eem',]
  p<-merge(a,b,by='key');d<-(p$estimate.x-p$truth.x)^2-(p$estimate.y-p$truth.y)^2
  report_close(if(length(d))mean(d)else NA_real_,paired$difference[i]);report_close(if(length(d)>1)sd(d)/sqrt(length(d))else NA_real_,paired$mcse[i])
  stopifnot(nrow(p)==paired$paired_n[i])
 }
 totals<-do.call(rbind,lapply(report_methods,function(mm){z<-r[r$method==mm,]
  data.frame(method=mm,attempted=nrow(z),point_n=sum(z$point_valid),ci_n=sum(z$interval_valid),
   point_failures=sum(!z$point_valid),interval_failures=sum(!z$interval_valid),rd_outside=sum(abs(z$estimate[z$point_valid])>1))}))
 coverage_lo<-min(pmax(0,s$coverage-1.96*s$coverage_mcse),na.rm=TRUE)
 limits<-list(rmse=range(pmax(.Machine$double.eps,s$rmse-1.96*s$rmse_mcse),s$rmse+1.96*s$rmse_mcse,na.rm=TRUE),
  bias=range(s$bias-1.96*s$bias_mcse,s$bias+1.96*s$bias_mcse,na.rm=TRUE),coverage=c(min(.88,floor(coverage_lo*20)/20),1))
 list(replicates=r,summary=s,paired=paired,totals=totals,limits=limits,
  settings=lapply(report_paths,function(p)read.csv(file.path(root,p,'settings.csv'))),checks=do.call(rbind,checks))
}
report_scenario_plot <- function(data,group) {
 z<-data$summary[data$summary$group==group,];z$method<-factor(z$method,levels=report_methods)
 z$sample_size<-factor(z$n,levels=c(500,1000),labels=c('n = 500','n = 1,000'));panels<-list()
 for(metric in c('rmse','bias','coverage')) {
  x<-z;x$value<-x[[metric]];x$lo<-x$value-1.96*x[[paste0(metric,'_mcse')]];x$hi<-x$value+1.96*x[[paste0(metric,'_mcse')]]
  if(metric=='rmse')x$lo<-pmax(.Machine$double.eps,x$lo)
  if(metric=='coverage'){x$lo<-pmax(0,x$lo);x$hi<-pmin(1,x$hi)}
  offset<-min(diff(sort(unique(x$level))))*.045
  p<-ggplot2::ggplot(x,ggplot2::aes(level,value,color=method,shape=method,linetype=method,group=method))+
   ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),width=offset,linewidth=.3,position=ggplot2::position_dodge(offset))+
   ggplot2::geom_line(linewidth=.5,position=ggplot2::position_dodge(offset))+
   ggplot2::geom_point(size=1.7,position=ggplot2::position_dodge(offset))+
   ggplot2::facet_wrap(~sample_size,nrow=1)+
   ggplot2::scale_color_manual(values=report_colors,labels=report_labels)+
   ggplot2::scale_shape_manual(values=c(original=15,direct_m_gam=16,drcmd=18,drcmd_eem=17),labels=report_labels)+
   ggplot2::scale_linetype_manual(values=c(original='dashed',direct_m_gam='solid',drcmd='dotted',drcmd_eem='dotdash'),labels=report_labels)+
   ggplot2::scale_x_continuous(breaks=sort(unique(x$level)))+
   ggplot2::labs(x=if(metric=='coverage')if(group==2)'Outcome slope b'else'Observation-probability scale s'else NULL,
    y=c(rmse='RMSE',bias='Bias',coverage='95% CI coverage')[[metric]],color=NULL,shape=NULL,linetype=NULL)+
   ggplot2::theme_classic(base_size=10)+ggplot2::theme(strip.background=ggplot2::element_blank(),strip.text=ggplot2::element_text(face='bold'),
    legend.position='bottom',panel.grid.major.y=ggplot2::element_line(color='grey90',linewidth=.2),axis.line=ggplot2::element_line(linewidth=.25))
  if(metric=='rmse')p<-p+ggplot2::scale_y_log10(limits=data$limits$rmse,breaks=c(.05,.1,.3,1,3,10,30),labels=c('0.05','0.1','0.3','1','3','10','30'))
  if(metric=='bias')p<-p+ggplot2::geom_hline(yintercept=0,linetype=2,color='grey55',linewidth=.3)+
   ggplot2::scale_y_continuous(trans=scales::pseudo_log_trans(sigma=.02),limits=data$limits$bias,breaks=c(-1,-.1,0,.1,1))
  if(metric=='coverage')p<-p+ggplot2::geom_hline(yintercept=.95,linetype=2,color='grey55',linewidth=.3)+
   ggplot2::scale_y_continuous(limits=data$limits$coverage,breaks=c(.90,.95,1))
  if(metric!='coverage')p<-p+ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
  if(metric!='rmse')p<-p+ggplot2::theme(strip.text=ggplot2::element_blank(),strip.background=ggplot2::element_blank())
  panels[[metric]]<-p
 }
 patchwork::wrap_plots(panels,ncol=1,guides='collect') & ggplot2::theme(legend.position='bottom',legend.text=ggplot2::element_text(size=8))
}
report_cell <- function(data,id,method,field) {
 z<-data$summary[data$summary$scenario==id&data$summary$method==method,field];stopifnot(length(z)==1);z
}
report_appendix_table <- function(data,group) {
 z<-data$summary[data$summary$group==group,];z<-z[order(z$n,z$level,match(z$method,report_methods)),]
 fmt<-function(x,se)sprintf('%.3f (%.3f)',x,se)
 tab<-data.frame(n=format(z$n,big.mark=',',trim=TRUE),Parameter=z$level,
  Estimator=unname(c(original='Ratio GAM',direct_m_gam='Direct GAM',drcmd='DRCMD',drcmd_eem='DRCMD EEM')[z$method]),
  'Bias (MCSE)'=fmt(z$bias,z$bias_mcse),'RMSE (MCSE)'=fmt(z$rmse,z$rmse_mcse),
  'Coverage (MCSE)'=fmt(z$coverage,z$coverage_mcse),'CI n'=z$ci_n,check.names=FALSE)
 names(tab)[2]<-if(group==2)'b'else's'
 knitr::kable(tab,format='latex',booktabs=TRUE,longtable=TRUE,row.names=FALSE,
  caption=paste(c('Observation-probability heterogeneity','Outcome slope','Quadratic outcome logit')[group],': numerical results.',sep=''))
}
