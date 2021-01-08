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
exp_s <- exp(-gamma*small_stepsize)
exp_b <- exp(-gamma*big_stepsize)

exp_s_gamma <- (1 - exp_s)/gamma
exp_b_gamma <- (1 - exp_b)/gamma

# Matrices
first_integral_var <- diag((1 - exp_b)/(2*gamma), dimension, dimension)
cov_const <- (1-exp_s) / (gamma*(1 + exp_s))
second_integral_var <- diag(((4*exp_s - exp_b + 2*gamma*small_stepsize - 3)/(2*(gamma^3))) -
                               (((1 - exp_s)^3) / (2*(gamma^3)*(1+exp_s))), dimension, dimension)

# Distribution of (x_0, v_0).
init <- function() {
  x <- t(fast_rmvnorm(1, mean = zero_vec, covariance = sigma_prior))
  v <- t(fast_rmvnorm(1, mean = zero_vec, covariance = u*identity_matrix))

  return(list(position = x, momentum = v, gradposition = gradlogtarget(x)))
}

chain_init <- function() {
  chain <- init()
  return(list(small_step_chain = chain, big_step_chain = chain))
}

# Small step of the Strang method
Strang_small_step <- function(chain_state, first_integral, second_integral) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradposition <- chain_state$gradposition
  
  v1 <- v0 + half_u*small_stepsize*gradposition
  v2 <- exp_s*v1 + sigma*first_integral
  
  x1 <- x0 + exp_s_gamma*v1 + sigma*second_integral

  gradposition <- gradlogtarget(x1)
  
  v3 <- v2 + half_u*small_stepsize*gradposition
  
  return(list(position = x1,
              momentum = v3,
              gradposition = gradposition))
}

# Big step of the Strang method
Strang_big_step <- function(chain_state, first_integral, second_integral) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradposition <- chain_state$gradposition
  
  v1 <- v0 + half_u*big_stepsize*gradposition
  v2 <- exp_b*v1 + sigma*first_integral
  
  x1 <- x0 + exp_b_gamma*v1 + sigma*second_integral
  
  gradposition <- gradlogtarget(x1)
  
  v3 <- v2 + half_u*big_stepsize*gradposition
  
  return(list(position = x1,
              momentum = v3,
              gradposition = gradposition))
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
  
  small_step_chain <- Strang_small_step(small_step_chain, first_integral_1, second_integral_1)
  small_step_chain <- Strang_small_step(small_step_chain, first_integral_2, second_integral_2)
  
  first_big_integral <- exp_s*first_integral_1 + first_integral_2
  
  second_big_integral <- second_integral_1 + exp_s_gamma * first_integral_1 + second_integral_2
  
  big_step_chain <- Strang_big_step(big_step_chain, first_big_integral, second_big_integral)
  
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