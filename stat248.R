### Load libraries and helper functions ###
source("helpers.R")

# Initialize constants for rows and columns
k = 3
j = 4

### Initialize Data Generating values of b ###
# true b's
bstrue = cbind(c(-2,0,1), c(3,1,0), c(0,-2,2), c(0,0,0), c(0,0,0))
b0 = bstrue[,1]
bs = bstrue[,2:5]

# Draw predictors
set.seed(100)
Xpredictors = mvtnorm::rmvnorm(n = 150, mean = c(0,0,0,0), sigma =diag(4) + ones(4))

# Get classifications
set.seed(131)
z = get_zs(b0, bs, Xpredictors)
C = apply(z, 1, function(x) which(rmultinom(1,1, get_prob(1:3, x)) == 1))
C

# Scatterplot with categoric color scale
plot(
  x = Xpredictors[,1],
  y = Xpredictors[,2],
  bg = C,
  cex = 1,
  pch=21
)


#### Full Logistic Regression ####
# Step 0: Set constants 
alpha_taustar = 1
beta_taustar = 1
epsilon = 0.05
alpha = 0.97

# Set seed 
set.seed(29)

# Draw from taustar
taustar = rexp(1, 1)
taus = rexp(j, taustar)
b0 = rnorm(k)
bs = matrix(rnorm(k*j, mean = 0, sd = 1/sqrt(rep( taus , k ))), nrow = k)
bfull = cbind(b0, bs)

# uncomment after first run to start second run circularly 
# bfull = matrix(btrace_mc[,100], nrow = k)
# bfull

Nsims = 100
btrace_mc = matrix(NA, nrow = 15, ncol = Nsims)
btrace_mc

#### Simulator takes around 80 min to run, Skip to loaded .csvs for results ####
for (iteration in 1:Nsims){
  print(paste("Iteration", iteration))
  for (step1rep in 1:10) {
    # Step 1a) 10 sims langevin
    p = draw_n() # Step 5
    
    btrace = matrix(NA, nrow = 15, ncol = 12)
    btrace[,1] = as.vector(bfull)
    btrace
    
    # Initialize p's
    ptrace = matrix(NA, nrow = 15, ncol = 12)
    ptrace[,1] = as.vector(p)
    ptrace
    
    
    results = run_mala_nsims(10, btrace, ptrace, taus)
    
    btrace = results$btrace
    ptrace = results$ptrace
    
    print("Mala done!")
    
    ### Step 1b) 25 Random-grid Metropolis updates for log(tau^*) using a proposal with w = 0.1
    
    posterior_log_target_taus <- function(x, taus, alpha, beta) dgamma(x, alpha + j, beta + sum(taus), log= T)
    posterior_log_specific <- function(x) posterior_log_target_taus(x, taus, alpha_taustar, beta_taustar)
    alpha_taustar = alpha_taustar + j 
    beta_taustar = beta_taustar + sum(taus)
    taustar <- rgm_final(25, posterior_log_specific, 0.1, taustar)
    
    
    ### Step 1c) Gibbs sampling updates for tau_j's
    draw_gibbs = function(alpha, beta, bs) {
      taus = c()
      betas = numeric(j)
      for(i in 1:j){
        taus = c(taus, rgamma(1, alpha + k/2, beta[i] + sum(bs[(1+i*k):(k+i*k)]^2)/2))
        betas[i] = beta[i] + sum(bs[(1+i*k):(k+i*k)]^2)/2
      }
      
      return(list(taus = taus, 
                  alpha = alpha + k/2, 
                  beta = betas))
    }
    result_gibbs = draw_gibbs( taustar, rep(1, j), btrace[,11])
    taus = result_gibbs$taus; taus 
    
    # current parameters 
    taus_alpha = result_gibbs$alpha; taus_alpha
    taus_beta = result_gibbs$beta; taus_beta
  }
  print("Step 1 completed!")
    
  ### Step 2) RGM for all the bjk simultaneously, 
  w = 0.01
  target_density = function(x) exp(log_target_density(matrix(x[4:15], nrow = 3), x[1:3]))
  results_rgm_multi <- rgm_final_multi(1, target_density, 0.01, btrace[,11])
  
  
  # Step 3) RGM update for log(tau*)
  w = 0.1
  # posterior_log_target_taus <- function(x, taus, alpha, beta) dgamma(x, alpha + j, beta + sum(taus), log= T)
  posterior_log_specific <- function(x) posterior_log_target_taus(x, taus, alpha_taustar, beta_taustar)
  
  taustar <- rgm_final(25, posterior_log_specific, 0.1, taustar)
  
  
  # Step 4) Gibbs sampling update for each of the taus
  result_gibbs = draw_gibbs( taustar, rep(1, j), results_rgm_multi)
  taus = result_gibbs$taus; taus ## CHANGE TO TAUS 
  
  # current parameters 
  taus_alpha = result_gibbs$alpha; taus_alpha
  taus_beta = result_gibbs$beta; taus_beta
  
  # Step 5) 
  p = rnorm(15)
  p
  btrace_mc[, iteration] = results_rgm_multi
  bfull = results_rgm_multi
  print(btrace_mc)
}

# B
Ytrace = read.csv('seed29_secondrun.csv')
Ytrace = as.matrix(Ytrace[,-1])
Xtrace = read.csv('seed29.csv')
Xtrace = as.matrix(Xtrace[,-1])


# Note faithful coupling after 85
for ( i in 1:100 ){
  if( sum(Ytrace[,i] == Xtrace[,i]) == 15){
    print(paste("Faithful coupling at step", i ))
    break
  }
}

# Plot all the individual components
for (i in 1:15) { 
  
  df = data.frame(t = rep(1:100, 2), val = c(Xtrace[i,],Ytrace[i,]), 
                  type = c(rep("X", 100), rep("Y", 100)))
  
  print(ggplot(df, aes(t, val, color = type)) + 
    geom_line() + 
    geom_vline( xintercept = 85) + 
    labs(title = paste("Component", i)) + 
    xlab("t") + 
    ylab("val") ) 
  Sys.sleep(2)
}


