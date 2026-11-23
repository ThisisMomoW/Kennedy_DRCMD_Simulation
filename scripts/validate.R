source('R/load.R');m<-check_study();start_observers()
# Compare numeric contents independently of score-column labels.
close_numeric<-function(x,y,tol=1e-10){
 x<-as.numeric(x);y<-as.numeric(y)
 stopifnot(length(x)==length(y),identical(is.na(x),is.na(y)))
 ok<-!is.na(x);stopifnot(all(abs(x[ok]-y[ok])<=tol*(1+abs(y[ok]))))
}
# Independent literal score calculation with both observed and missing exposures.
d<-data.frame(X1=0,X2=0,Y=c(0,1,0,1,1),R=c(1,1,0,0,1),A_obs=c(0,1,NA,NA,0))
k<-c(.2,.8,.3,.4,.5);l<-c(.1,.8,.2,.7,.4);g<-c(.4,.6,.3,.7,.5);b0<-(1-g)*.3;b1<-g*.6
reference<-function(d,k,l,g,b0,b1){
  ans<-matrix(0,nrow(d),2)
  for(i in seq_len(nrow(d)))for(a in 0:1){
    ga<-if(a==1)g[i]else 1-g[i];la<-if(a==1)l[i]else 1-l[i]
    ma<-if(a==1)b1[i]/ga else b0[i]/ga
    correction<-if(d$R[i]==1)((d$A_obs[i]==a)-la)/k[i]else 0
    ans[i,a+1]<-ma+(d$Y[i]-ma)/ga*(la+correction)
  };ans
}
x<-kennedy_score(d,k,l,g,b0,b1);close_numeric(x,reference(d,k,l,g,b0,b1))
z<-kennedy_finish(x);close_numeric(z$estimate,mean(x[,2]-x[,1]));close_numeric(z$se,sqrt(var(x[,2]-x[,1])/nrow(d)))
o<-c(5,2,4,1,3);close_numeric(kennedy_score(d[o,],k[o],l[o],g[o],b0[o],b1[o]),x[o,])
d$latent_A<-0;close_numeric(kennedy_score(d,k,l,g,b0,b1),x);d$latent_A<-1;close_numeric(kennedy_score(d,k,l,g,b0,b1),x)
d$R<-1;d$A_obs<-c(0,1,0,1,0)
y<-kennedy_score(d,rep(1,5),l,g,b0,b1)
close_numeric(y[,1],.3+(d$A_obs==0)*(d$Y-.3)/(1-g));close_numeric(y[,2],.6+(d$A_obs==1)*(d$Y-.6)/g)
# Numerical fidelity on an independently saved reference fixture.
f<-readRDS('reference/verification_fixture.rds');study_context$method<-'drcmd';study_audit_reset()
e<-get('.gam_audit',.GlobalEnv);e$active<-TRUE;set_stream(f$drcmd_seed);dr<-drcmd_point(f$observed);e$active<-FALSE
stopifnot(identical(dr$audit$folds,f$folds),dr$audit$counts['drcmd']==1,
          dr$audit$counts['est_varphi_eem']==0,dr$audit$counts['est_psi_tml']==0)
for(z in dr$audit$sl){stopifnot(identical(z$library,'SL.gam'));if(z$stage=='est_varphi')stopifnot(z$family=='gaussian')}
set_stream(f$kennedy_seed);study_context$method<-'kennedy';study_context$fold<-0L
kr<-capture_fit(kennedy_fit(f$observed,f$folds));stopifnot(!is.null(kr$value))
kd<-kennedy_direct_fit(f$observed,f$folds,kr);stopifnot(!is.null(kd$result))
checks<-list(drcmd=dr,original=kr$value$result,direct_m_gam=kd$result)
for(mm in names(checks))close_numeric(unlist(checks[[mm]][c('estimate','se','lower','upper')]),as.numeric(f$expected[f$expected$method==mm,c('estimate','se','lower','upper')]))
# Confirm the direct model's fractional weights, formula and prediction scale independently.
ff<-f$folds[[1]];p<-kr$value$fold_details[[1]];tr<-f$observed[ff$train,];te<-f$observed[ff$test,]
for(a in 0:1){
  w<-if(a==1)p$generated_training$lambda1 else 1-p$generated_training$lambda1
  tr$case_weight<-w
  manual<-gam::gam(Y~s(X1,2)+s(X2,2),data=tr,family=quasibinomial('logit'),weights=case_weight,
                   control=gam::gam.control(maxit=50,bf.maxit=50))
  actual<-kd$fold_details[[1]]$fits[[a+1]]
  close_numeric(as.numeric(predict(manual,te,type='response')),actual$value$prediction)
  stopifnot(identical(actual$gam_fit_audit$input_weights,as.numeric(w)),
            identical(actual$value$prior_weights,as.numeric(w)),any(w>0&w<1),sum(tr$R==0)>0)
}
alt<-te;alt$Y<-1-alt$Y;alt$A_obs<-1
z<-dmg_fit_fold(tr,alt,p$generated_training$lambda1)
close_numeric(z$m0,kd$fold_details[[1]]$m0);close_numeric(z$m1,kd$fold_details[[1]]$m1)
for(piece in kr$value$fold_details)for(stage in piece$fit_audit)stopifnot(!any(stage$fit_ids%in%piece$test_ids))
# Verify Bayes' rule and identities for binary Y over an independent quadrature grid.
xx<-quad_grid(40)
for(i in seq_len(nrow(m$settings))){
  s<-m$settings[i,];p<-dgp_function(s$group)(xx,s);py<-p$g*p$m1+(1-p$g)*p$m0
  lam1<-p$g*p$m1/py;lam0<-p$g*(1-p$m1)/(1-py)
  close_numeric(py*lam1+(1-py)*lam0,p$g);close_numeric(py*lam1,p$g*p$m1)
  close_numeric(sum(xx$weight*(p$m1-p$m0)),s$truth)
}
writeLines(c('PASS: reference estimates, SEs and intervals reproduced for all three methods.',
 'PASS: population calibration, binary Bayes identities, marginal RD, and quadrature accuracy.',
 'PASS: native exported drcmd call, Gaussian augmentation, native folds, no EEM or TMLE.',
 'PASS: score identity, complete-data reduction, missing-value safety, permutation and paired SE.',
 'PASS: fold isolation; fractional weights at gam.fit; manual gam formula and predictions.'),'validation/VERIFIED')
write.csv(f$expected,'validation/reference_fixture_results.csv',row.names=FALSE)
cat(readLines('validation/VERIFIED'),sep='\n')
