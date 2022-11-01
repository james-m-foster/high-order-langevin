source("GermanCredit.R")

set.seed(1)

# Prior covariance for regression coefficients
sigma2 <- 10
sigma_prior <- diag(sigma2, dimension, dimension)
identity_matrix <- diag(1, dimension, dimension)
zero_vec <- matrix(0, nrow = dimension, ncol = 1)
logistic_setting <- logisticregression_precomputation(Y, X, zero_vec, sigma_prior)

# ULD parameters
gamma <- 2
u <- 1
sigma <- sqrt(2*gamma*u)

half_u <- 0.5*u

# Step size
big_stepsize <- 0.01
small_stepsize <- 0.5*big_stepsize

# Precomputed constants
exp_hs <- exp(-0.5*gamma*small_stepsize)
exp_s <- exp(-gamma*small_stepsize)
exp_b <- exp(-gamma*big_stepsize)

exp_hs_gamma <- (1 - exp_hs)/gamma
exp_s_gamma <- (1 - exp_s)/gamma
exp_b_gamma <- (1 - exp_b)/gamma

exp_s_gamma2 <- (exp_s + gamma * small_stepsize - 1)/(gamma^2)
exp_b_gamma2 <- (exp_b + gamma * big_stepsize - 1)/(gamma^2)

# Matrices
first_integral_var <- diag((1 - exp_s)/(2*gamma), dimension, dimension)
cov_const <- (1-exp_hs) / (gamma*(1 + exp_hs))
second_integral_var <- diag(((4*exp_hs - exp_s + gamma*small_stepsize - 3)/(2*(gamma^3))) -
                              (((1 - exp_hs)^3) / (2*(gamma^3)*(1+exp_hs))), dimension, dimension)

# Distribution of (x_0, v_0).
init <- function() {
  x <- t(fast_rmvnorm(1, mean = zero_vec, covariance = sigma_prior))
  v <- t(fast_rmvnorm(1, mean = zero_vec, covariance = u*identity_matrix))
  
  return(list(position = x, momentum = v))
}

chain_init <- function() {
  chain <- init()
  return(list(small_step_chain = chain, big_step_chain = chain))
}

# Small step of the UBU method (Sanz-Serna and Zygalakis (JMLR, 2021))
UBU_small_step <- function(chain_state, little_integral, first_integral, second_integral) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradiant <- gradlogtarget(x0 + exp_hs_gamma*v0 + sigma*little_integral)
  
  v1 <- exp_s*v0 + u*exp_s_gamma*gradiant + sigma*first_integral
  x1 <- x0 + exp_s_gamma*v0 + u*exp_s_gamma2*gradiant + sigma*second_integral
  
  return(list(position = x1,
              momentum = v1))
}

# Big step of the UBU method (Sanz-Serna and Zygalakis (JMLR, 2021))
UBU_big_step <- function(chain_state, little_integral, first_integral, second_integral) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradiant <- gradlogtarget(x0 + exp_s_gamma*v0 + sigma*little_integral)
  
  v1 <- exp_b*v0 + u*exp_b_gamma*gradiant + sigma*first_integral
  x1 <- x0 + exp_b_gamma*v0 + u*exp_b_gamma2*gradiant + sigma*second_integral
  
  return(list(position = x1,
              momentum = v1))
}

# Kernel for coupled ULD chain
coupled_kernel <- function(chains) {
  small_step_chain <- chains$small_step_chain
  big_step_chain <- chains$big_step_chain
  
    
  first_integral_1 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = first_integral_var))
  second_integral_1 <- t(fast_rmvnorm(1, mean = cov_const*first_integral_1,
                                      covariance = second_integral_var))
  
  first_integral_2 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = first_integral_var))
  second_integral_2 <- t(fast_rmvnorm(1, mean = cov_const*first_integral_2,
                                      covariance = second_integral_var))
  
  first_integral_3 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = first_integral_var))
  second_integral_3 <- t(fast_rmvnorm(1, mean = cov_const*first_integral_3,
                                      covariance = second_integral_var))
  
  first_integral_4 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = first_integral_var))
  second_integral_4 <- t(fast_rmvnorm(1, mean = cov_const*first_integral_4,
                                      covariance = second_integral_var))
  
  
  first_integral_5 <- exp_hs*first_integral_1 + first_integral_2
  
  second_integral_6 <- second_integral_1 + exp_hs_gamma * first_integral_1 + second_integral_2
  
  first_integral_7 <- exp_hs*first_integral_3 + first_integral_4
  
  second_integral_8 <- second_integral_3 + exp_hs_gamma * first_integral_3 + second_integral_4
  
  
  small_step_chain <- UBU_small_step(small_step_chain, second_integral_1 , first_integral_5, second_integral_6)
  small_step_chain <- UBU_small_step(small_step_chain, second_integral_3 , first_integral_7, second_integral_8)
  
  first_big_integral <- exp_s*first_integral_5 + first_integral_7
  
  second_big_integral <- second_integral_6 + exp_s_gamma * first_integral_5 + second_integral_8
  
  big_step_chain <- UBU_big_step(big_step_chain, second_integral_6, first_big_integral, second_big_integral)
  
  return(list(small_step_chain = small_step_chain,
              big_step_chain = big_step_chain))
}

# Run numerical experiment
nsamples <- 100
end_time <- 1000
no_of_steps <- end_time / big_stepsize
error <- 0

for (i in 1:nsamples) {
  chain <- chain_init()
  
  for (j in 1:no_of_steps) {
    chain <- coupled_kernel(chain)
  }
  
  error <- error + sum((chain$small_step_chain$position - chain$big_step_chain$position)^2)
}

error <- error / nsamples

print(sqrt(error))
