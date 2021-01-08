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
exp_s <- exp(-0.5*gamma*small_stepsize)
exp_b <- exp(-gamma*small_stepsize)

# Matrices
G_small_var <- diag(((1 - exp_b)/(2*gamma)), dimension, dimension)

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

# Small step of the OBABO method
OBABO_small_step <- function(chain_state, G_small_increment_1, G_small_increment_2) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradposition <- chain_state$gradposition
  
  v1 <- exp_s*v0 + sigma * G_small_increment_1 + half_u*small_stepsize*gradposition
  x1 <- x0 + v1*small_stepsize
  
  gradposition <- gradlogtarget(x1)
  
  v2 <- exp_s*(v1 + half_u*small_stepsize*gradposition) + sigma * G_small_increment_2 
  
  return(list(position = x1,
              momentum = v2,
              gradposition = gradposition))
}

# Big step of the OBABO method
OBABO_big_step <- function(chain_state, G_big_increment_1, G_big_increment_2) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  gradposition <- chain_state$gradposition
  
  v1 <- exp_b*v0 + sigma * G_big_increment_1 + half_u*big_stepsize*gradposition
  x1 <- x0 + v1*big_stepsize
  
  gradposition <- gradlogtarget(x1)
  
  v2 <- exp_b*(v1 + half_u*big_stepsize*gradposition) + sigma * G_big_increment_2

  
  return(list(position = x1,
              momentum = v2,
              gradposition = gradposition))
}

# Kernel for coupled ULD chain
coupled_kernel <- function(chains) {
  small_step_chain <- chains$small_step_chain
  big_step_chain <- chains$big_step_chain
  
  G_small_increment_1 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = G_small_var))
  G_small_increment_2 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = G_small_var))
    
  G_small_increment_3 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = G_small_var))
  G_small_increment_4 <- t(fast_rmvnorm(1, mean = zero_vec, covariance = G_small_var))
  
  small_step_chain <- OBABO_small_step(small_step_chain, G_small_increment_1, G_small_increment_2)
  small_step_chain <- OBABO_small_step(small_step_chain, G_small_increment_3, G_small_increment_4)
  
  G_big_increment_1 <- exp_s*G_small_increment_1 + G_small_increment_2
  
  G_big_increment_2 <- exp_s*G_small_increment_3 + G_small_increment_4
  
  big_step_chain <- OBABO_big_step(big_step_chain, G_big_increment_1, G_big_increment_2)
  
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
