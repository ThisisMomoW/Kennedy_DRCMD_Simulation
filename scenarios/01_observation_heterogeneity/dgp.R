# X1,X2 independent U(-1,1); only the observation scale s varies.
dgp_parameters <- function(x,setting) {
  list(g=plogis(.4*x$X1-.4*x$X2),
       m0=plogis(-.5+.5*x$X1+.5*x$X2),
       m1=plogis(.5*x$X1+.5*x$X2),scale=setting$level)
}
