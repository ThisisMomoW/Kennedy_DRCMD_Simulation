dgp_function<-function(group) {
  e<-new.env(parent=globalenv());sys.source(file.path(scenario_paths[group],'dgp.R'),envir=e);e$dgp_parameters
}
quad_grid <- function(N=80L) {
  j<-1:(N-1);J<-matrix(0,N,N);J[cbind(j,j+1)]<-J[cbind(j+1,j)]<-j/sqrt(4*j*j-1)
  e<-eigen(J,symmetric=TRUE);z<-expand.grid(X1=e$values,X2=e$values)
  z$weight<-as.vector(outer(e$vectors[1,]^2,e$vectors[1,]^2));z
}
settings_grid <- function() {
  z<-do.call(rbind,lapply(1:3,function(g){
    x<-expand.grid(n=c(500L,1000L),level=switch(g,c(1,3,5),c(.5,1.5,3),c(1,5)))
    x$group<-g;x$delta<-if(g==3)2 else 0;x
  }))
  z$scenario<-c(sprintf('e%02d',1:12),sprintf('nl%02d',1:4));z$q<-.8;z$attempted<-200L;z
}
population_settings<-function(N=80L) {
  grid<-quad_grid(N);ss<-settings_grid()
  do.call(rbind,lapply(seq_len(nrow(ss)),function(i){
    s<-ss[i,];f<-dgp_function(s$group)(grid,s);p<-f$g*f$m1+(1-f$g)*f$m0
    # Sum exactly over A,Y, using P(Y=1|X)=p after marginalizing A.
    rate<-function(c)sum(grid$weight*((1-p)*plogis(c+f$scale*.4*grid$X1)+p*plogis(c+f$scale*(.4*grid$X1+.2))))
    cc<-uniroot(function(c)rate(c)-(1-s$q),c(-20,20),tol=1e-13)$root
    data.frame(s,c=cc,truth=sum(grid$weight*(f$m1-f$m0)),risk0=sum(grid$weight*f$m0),
               risk1=sum(grid$weight*f$m1),observation_rate=rate(cc),integration_N=N)
  }))
}
generate_dataset<-function(setting) {
  n<-setting$n;x<-data.frame(X1=runif(n,-1,1),X2=runif(n,-1,1));f<-dgp_function(setting$group)(x,setting)
  a<-rbinom(n,1,f$g);x$Y<-rbinom(n,1,ifelse(a==1,f$m1,f$m0))
  kappa<-plogis(setting$c+f$scale*(.4*x$X1+.2*x$Y));x$R<-rbinom(n,1,kappa)
  x$A_obs<-ifelse(x$R==1,a,NA_real_)
  p<-f$g*f$m1+(1-f$g)*f$m0
  lambda<-f$g*ifelse(x$Y==1,f$m1,1-f$m1)/ifelse(x$Y==1,p,1-p)
  list(observed=x,evaluation=list(latent_A=a,g=f$g,m0=f$m0,m1=f$m1,kappa=kappa,lambda=lambda))
}
