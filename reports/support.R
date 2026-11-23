# Saved-result checks and graphics only; no estimator is called during rendering.
report_methods<-c('original','direct_m_gam','drcmd')
report_labels<-c(original='Kennedy ratio GAM',direct_m_gam='Kennedy direct GAM',drcmd='DR-CMD')
report_colors<-c(original='#B8411A',direct_m_gam='#00867D',drcmd='#2E5378')
report_paths<-c('scenarios/01_observation_heterogeneity','scenarios/02_outcome_slope','scenarios/03_nonlinear_outcome')
read_report_data<-function(root){
 r<-read.csv(file.path(root,'results/replicates.csv'));s<-read.csv(file.path(root,'results/summary.csv'))
 stopifnot(nrow(r)==9600,nrow(s)==48,length(unique(r$key))==3200,all(s$attempted==200),
           setequal(unique(r$method),report_methods),!anyDuplicated(r[c('key','method')]))
 checks<-list()
 for(i in seq_len(nrow(s))){
  z<-r[r$scenario==s$scenario[i]&r$method==s$method[i],];e<-z$estimate[is.finite(z$estimate)]-z$truth[is.finite(z$estimate)]
  ci<-z[z$interval_valid,];covered<-ci$lower<=ci$truth&ci$truth<=ci$upper
  rr<-sqrt(mean(e^2));cc<-mean(covered)
  vals<-c(bias=mean(e),bias_mcse=sd(e)/sqrt(length(e)),rmse=rr,rmse_mcse=sd(e^2)/(2*rr*sqrt(length(e))),
    coverage=cc,coverage_mcse=sqrt(cc*(1-cc)/length(covered)))
  err<-max(abs(vals-as.numeric(s[i,names(vals)])))
  stopifnot(err<1e-10*(1+max(abs(vals))),nrow(z)==200,length(e)==s$point_n[i],nrow(ci)==s$ci_n[i])
  checks[[i]]<-data.frame(scenario=s$scenario[i],method=s$method[i],max_difference=err)
 }
 for(z in split(r,r$key))stopifnot(length(unique(z$data_hash))==1,length(unique(z$fold_hash))==1)
 for(g in 1:3){loc<-read.csv(file.path(root,report_paths[g],'results/summary.csv'));a<-s[s$group==g,]
  b<-loc[match(paste(a$scenario,a$method),paste(loc$scenario,loc$method)),]
  stopifnot(isTRUE(all.equal(a$rmse,b$rmse)),isTRUE(all.equal(a$coverage,b$coverage)))
 }
 paired<-read.csv(file.path(root,'results/paired_mse.csv'))
 for(i in seq_len(nrow(paired))){
  z<-r[r$scenario==paired$scenario[i]&r$point_valid,];a<-z[z$method=='drcmd',];b<-z[z$method==paired$comparator[i],]
  p<-merge(a,b,by='key');d<-(p$estimate.x-p$truth.x)^2-(p$estimate.y-p$truth.y)^2
  stopifnot(nrow(p)==paired$paired_n[i],abs(mean(d)-paired$difference[i])<1e-9*(1+abs(mean(d))),
    abs(sd(d)/sqrt(length(d))-paired$mcse[i])<1e-9*(1+paired$mcse[i]))
 }
 totals<-do.call(rbind,lapply(report_methods,function(mm){z<-r[r$method==mm,]
  data.frame(method=mm,attempted=nrow(z),point_n=sum(z$point_valid),ci_n=sum(z$interval_valid),
   point_failures=sum(!z$point_valid),interval_failures=sum(!z$interval_valid),rd_outside=sum(abs(z$estimate[z$point_valid])>1))}))
 coverage_lo<-min(pmax(0,s$coverage-1.96*s$coverage_mcse))
 limits<-list(rmse=range(pmax(.Machine$double.eps,s$rmse-1.96*s$rmse_mcse),s$rmse+1.96*s$rmse_mcse),
  bias=range(s$bias-1.96*s$bias_mcse,s$bias+1.96*s$bias_mcse),coverage=c(min(.88,floor(coverage_lo*20)/20),1))
 list(replicates=r,summary=s,paired=paired,totals=totals,limits=limits,
  settings=lapply(report_paths,function(p)read.csv(file.path(root,p,'settings.csv'))),
  checks=do.call(rbind,checks),reference=read.csv(file.path(root,'results/reference_comparison.csv')))
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
   ggplot2::scale_shape_manual(values=c(original=15,direct_m_gam=16,drcmd=18),labels=report_labels)+
   ggplot2::scale_linetype_manual(values=c(original='dashed',direct_m_gam='solid',drcmd='dotted'),labels=report_labels)+
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
 patchwork::wrap_plots(panels,ncol=1,guides='collect') & ggplot2::theme(legend.position='bottom')
}
report_cell <- function(data,id,method,field) {
 z<-data$summary[data$summary$scenario==id&data$summary$method==method,field];stopifnot(length(z)==1);z
}
report_appendix_table <- function(data,group) {
 z<-data$summary[data$summary$group==group,];z<-z[order(z$n,z$level,match(z$method,report_methods)),]
 fmt<-function(x,se)sprintf('%.3f (%.3f)',x,se)
 tab<-data.frame(n=format(z$n,big.mark=',',trim=TRUE),Parameter=z$level,
  Estimator=unname(c(original='Ratio GAM',direct_m_gam='Direct GAM',drcmd='DR-CMD')[z$method]),
  'Bias (MCSE)'=fmt(z$bias,z$bias_mcse),'RMSE (MCSE)'=fmt(z$rmse,z$rmse_mcse),
  'Coverage (MCSE)'=fmt(z$coverage,z$coverage_mcse),'CI n'=z$ci_n,check.names=FALSE)
 names(tab)[2]<-if(group==2)'b'else's'
 knitr::kable(tab,format='latex',booktabs=TRUE,longtable=TRUE,row.names=FALSE,
  caption=paste(c('Observation-probability heterogeneity','Outcome slope','Quadratic outcome logit')[group],': numerical results.',sep=''))
}
