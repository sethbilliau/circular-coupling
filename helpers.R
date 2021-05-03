### Helper Functions 
library(mvtnorm)
library(matlab)
library(rootSolve)
library(paletteer)


#### Data Generating ####

# Function to get z given predictors
get_zs = function(b0, bs, predictors = Xpredictors) {
  z = matrix(NA, nrow = 150, ncol = 3)
  for (k in 1:3){
    for (i in 1:150){
      z[i,k]=  b0[k] + sum(bs[k,] * predictors[i,])
    }
  }
  return(z)
}

# Function to get probability given Z's
get_prob = function(i, vector) exp(vector[i])/sum(exp(vector))

# plot
colors <- paletteer_c(package = "ggthemes", palette = "Green-Blue-White", n = 3)
colors

#### Functions for Logistic Regression ####
# Get log likelihood
logLik = function(z, Cs=C) {
  loglikelihood = sum(apply(cbind(C, z), 1, function(x) log(get_prob(x[1], x[2:4]))))
  return(loglikelihood)
}

# Log prior
logPrior= function(taus, bs, b0) { 
  logprior=0
  for (i in 1:length(bs)){
    logprior = logprior + dnorm(bs[i],mean =0, sd= 1 / sqrt(taus[ceil(i/3)]), log= T)
  }
  for (i in 1:length(b0)){ 
    logprior = logprior + dnorm(b0[i],mean =0, sd= 1, log= T)
  }
  return(logprior)
} 

# Get log target density
log_target_density = function(bs, b0, tau=taus, predictors = Xpredictors, Cs = C) { 
  z = get_zs(b0, bs, predictors)
  result = logLik(z, C) + logPrior(tau, bs, b0)
  return(result)
}



#### MALA ####
# Run Mala once
run_mala = function(bfull, p, taus){
  # Replace p
  p = update_p(p)
  
  # Get p tick
  Efun = function(x) -log_target_density( x[,2:5],x[,1], taus)
  get_p_tick = function(p, x) p - matrix((epsilon/2)*gradient(Efun, x), nrow=3)
  ptick = get_p_tick(p, bfull)
  
  # Set xstar
  bfullstar = get_bfullstar(bfull, ptick)
  
  # Set pstar 
  get_pstar = function(bfullstar_old, ptick_old) ptick_old - matrix((epsilon/2)*gradient(Efun, bfullstar_old), nrow=3)
  pstar = get_pstar(bfullstar, ptick)
  
  # Accept/Reject
  Hfun = function(b, p) Efun(b) + sum(as.vector(p)^2)/2
  U = runif(1)
  acceptprob = min(1, exp(Hfun(bfull, p) - Hfun(bfullstar, pstar)))
  if (U < acceptprob) { 
    return(list(x= as.vector(bfullstar), p=pstar, accept=1))
  } else {
    return(list(x= as.vector(bfull), p=-p, accept=0))
  }
}

# Run MALA n times
run_mala_nsims = function(nsims, btrace, ptrace, taus) {
  acceptcounter = 0
  for (i in 1:nsims){
    bmat = matrix(btrace[,i], nrow=3)
    pmat = matrix(ptrace[,i], nrow=3)
    result= run_mala(bmat, pmat, taus)
    btrace[,i+1] = result$x
    ptrace[,i+1] = result$p
    acceptcounter = acceptcounter + result$accept
    if(i %% 100 == 0){
      print(i)
    }
  }
  
  return(list(btrace = btrace, 
              ptrace = ptrace, 
              accept = acceptcounter/nsims))
}




#### Random Grid code ####
f <- function(w, x, u) {
  return(2*w * ((u[2] - 1/2) + round(x / (2*w) - (u[2] - 1/2))))
}



phi <- function(pi, w, x, u) {
  # print( pi(x))
  if (u[1] < (pi(f(w, x, u)) / pi(x)) ){
    return(f(w, x, u))
  }
  return(x)
}

rgm_history <- function(n, pi, w, x0, unifs = NULL) {
  if (is.null(unifs)) {
    unifs <- replicate(n, runif(2))
  }
  
  states <- numeric(n + 1)
  states[1] <- x0
  for (i in 1:n) {
    states[i+1] <- phi(pi, w, states[i], unifs[, i])
  }
  
  return(list(states, unifs))
}

rgm_final <- function(n, pi, w, x0, unifs = NULL) {
  return(rgm_history(n, pi, w, x0, unifs)[[1]][n + 1])
}

f_multi <- function(w, x, u) {
  uis <- tail(u, length(x))
  return(2*w * ((uis - 1/2) + round(x / (2*w) - (uis - 1/2))))
}

phi_multi <- function(pi, w, x, u) {
  if (u[1] < pi(f(w, x, u)) / pi(x)) {
    return(f_multi(w, x, u))
  }
  return(x)
}

rgm_history_multi <- function(n, pi, w, x0, unifs = NULL) {
  if (is.null(unifs)) {
    unifs <- replicate(n, runif(length(x0) + 1))
  }
  
  states <- data.frame(replicate(n + 1, numeric(length(x0))))
  states[,1] = x0
  for (i in 1:n) {
    states[, i+1] <- phi_multi(pi, w, states[, i], unifs[, i])
  }
  
  return(list(states, unifs))
}

rgm_final_multi <- function(n, pi, w, x0, unifs = NULL) {
  return(rgm_history_multi(n, pi, w, x0, unifs)[[1]][, n + 1])
}


