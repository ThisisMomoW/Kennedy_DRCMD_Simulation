# X2 remains a supplied covariate although it has no generating effect.
dgp_parameters <- function(x,setting) {
  list(g=plogis(.4*x$X1),m0=plogis(-.5+setting$level*x$X1),
       m1=plogis(setting$level*x$X1),scale=1)
}
