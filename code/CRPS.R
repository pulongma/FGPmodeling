CRPS <- function(x, mu, sig){
  ### Input Arguments:
  # x: a vector of true values to be predicted
  # mu: a vector of predictive mean
  # sig: a vector of predictive standard deviation
  
  xo = (x-mu)/sig
  crps = sig*(xo*(2*pnorm(xo)-1) + 2*dnorm(xo) - 1/sqrt(pi))
  
  return(crps)
}
