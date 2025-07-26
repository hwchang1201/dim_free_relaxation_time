rm(list=ls())
# XtX = {
#   {n,   c*n, c*n, c*n, c*n, c*n},
#   {c*n,   n, c*n, c*n, c*n, c*n},
#   {c*n, c*n,   n, c*n, c*n, c*n},
#   {c*n, c*n, c*n,   n, c*n, c*n},
#   {c*n, c*n, c*n, c*n,   n, c*n},
#   {c*n, c*n, c*n, c*n, c*n,   n},
# };

# Xty = {b*n, b*c*n, b*c*n, b*c*n, b*c*n, b*c*n};
# yty = (b^2*n  + n);

n = 100; p = 6; c = 0.5; b = 5; 
kappa = 2; g = p^2
ell = p; L = p^3

XtX = matrix(c*n, nrow = p, ncol = p); diag(XtX) = n
Xty = rep(b*c*n, p); Xty[1] = b*n
yty = b^2*n + n

clip <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

logposterior <- function(model){
  size = sum(model)
  idx = which(model == 1)
  logposterior0 = -(n/2) * log(1+g)
  if(size == 0){
    return(0)
  }
  RSS = as.numeric(yty - crossprod(Xty[idx], solve(XtX[idx, idx], Xty[idx])))/yty
  logpost = - size* (kappa*log(p) + log(1+g)/2 ) -(n/2) * log(1+g * RSS)
  
  logpost - logposterior0
}

# These are the cases that we only consider.
logposterior(c(0,0,0,0,0,0))
logposterior(c(1,0,0,0,0,0))
logposterior(c(0,1,0,0,0,0))
logposterior(c(1,1,0,0,0,0))
logposterior(c(0,1,1,0,0,0))
logposterior(c(1,1,1,0,0,0))
log(ell)
log(L)

proposal_probability <- function(model1, model2){
  if(sum(abs(model1 - model2)) != 1){
    return(NA)
  }
  clipped_value = numeric(p)
  for (i in 1:p){
    model3 = model1
    model3[i] = 1 - model3[i]
    clipped_value[i] = clip(logposterior(model3) - logposterior(model1), log(ell), log(L))
    if(sum(abs(model3 - model2)) == 0){
      j = i 
    }
  }
  vec = exp(clipped_value - max(clipped_value))
  vec = vec/sum(vec)
  return(vec[j])
}

transition_probability <- function(model1, model2){
  if(sum(abs(model1 - model2)) != 1){
    return(NA)
  }
  prop12 = proposal_probability(model1, model2)
  prop21 = proposal_probability(model2, model1)
  
  prop12 * min(1, exp(logposterior(model2) - logposterior(model1)) * prop21 /prop12)
}

# empty -> {1} : prop = 1/p
transition_probability(c(0,0,0,0,0,0), c(1,0,0,0,0,0))
# empty -> {j} : prop = 1/p, for j = 2,3,4,5,6
transition_probability(c(0,0,0,0,0,0), c(0,1,0,0,0,0))
# {j} -> {1,j} : prop = p^3/(p^3+p^2), for j = 2,3,4,5,6
transition_probability(c(0,1,0,0,0,0), c(1,1,0,0,0,0))
# {1,j} -> {1} : prop = p^3/(p^3+p^2), for j = 2,3,4,5,6
transition_probability(c(1,1,0,0,0,0), c(1,0,0,0,0,0))


