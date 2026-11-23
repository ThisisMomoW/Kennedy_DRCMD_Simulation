# The quadratic term is in the outcome generator, not a supplied model feature.
dgp_parameters <- function(x,setting) {
  eta<-.5*x$X1+.5*x$X2+2*(x$X1^2-1/3)
  list(g=plogis(.4*x$X1-.4*x$X2),m0=plogis(-.5+eta),m1=plogis(eta),scale=setting$level)
}
