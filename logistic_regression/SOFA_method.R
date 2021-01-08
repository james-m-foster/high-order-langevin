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

# Step size
big_stepsize <- 0.01
big_stepsize_squared <- big_stepsize ** 2
no_of_small_steps <- 2
small_stepsize <- big_stepsize / no_of_small_steps
small_stepsize_squared <- small_stepsize ** 2
small_stepsize_halfed <- 0.5 * small_stepsize

# Precomputed constants
third <- 1/3

splitting_const <- ((2 ** third) - 1)/(2*(2 - (2 ** third)))

exp_1_s <- exp(-gamma*(0.5 + splitting_const)*small_stepsize)
exp_2_s <- exp(gamma*splitting_const*small_stepsize)
exp_integral_1_s <- (1 - exp_1_s)/(gamma*small_stepsize)
exp_integral_2_s <- (1 - exp_2_s)/(gamma*small_stepsize)

exp_1_b <- exp(-gamma*(0.5 + splitting_const)*big_stepsize)
exp_2_b <- exp(gamma*splitting_const*big_stepsize)
exp_integral_1_b <- (1 - exp_1_b)/(gamma*big_stepsize)
exp_integral_2_b <- (1 - exp_2_b)/(gamma*big_stepsize)

# Matrices
W_small_var <- diag(small_stepsize, dimension, dimension)
H_small_var <- diag(small_stepsize / 12, dimension, dimension)
K_small_var <- diag(small_stepsize / 720, dimension, dimension)

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

# Small step of the proposed method
ULD_small_step <- function(chain_state, W_small_increment, H_small_increment, K_small_increment) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum + sigma*(H_small_increment + 6*K_small_increment)
  
  gradposition <- chain_state$gradposition
  
  diagonal_term <- sigma*(W_small_increment - 12*K_small_increment)
  
  v1 <- exp_1_s*v0 + (u * small_stepsize * gradposition + diagonal_term)*exp_integral_1_s
  x1 <- x0 + v1*(1+2*splitting_const)*small_stepsize
  
  v2 <- exp_2_s*v1 + (u * small_stepsize * gradlogtarget(x1) + diagonal_term)*exp_integral_2_s
  x2 <- x1 - v2*(1+4*splitting_const)*small_stepsize
  
  v3 <- exp_2_s*v2 + (u * small_stepsize * gradlogtarget(x2) + diagonal_term)*exp_integral_2_s
  x3 <- x2 + v3*(1+2*splitting_const)*small_stepsize
  
  gradposition <- gradlogtarget(x3)
  
  v4 <- exp_1_s*v3 + (u * small_stepsize * gradposition + diagonal_term)*exp_integral_1_s

  return(list(position = x3,
              momentum = v4 - sigma*(H_small_increment - 6*K_small_increment),
              gradposition = gradposition))
}

# Big step of the proposed method
ULD_big_step <- function(chain_state, W_big_increment, H_big_increment, K_big_increment) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum + sigma*(H_big_increment + 6*K_big_increment)
  
  gradposition <- chain_state$gradposition
  
  diagonal_term <- sigma*(W_big_increment - 12*K_big_increment)
  
  v1 <- exp_1_b*v0 + (u * big_stepsize * gradposition + diagonal_term)*exp_integral_1_b
  x1 <- x0 + v1*(1+2*splitting_const)*big_stepsize
  
  v2 <- exp_2_b*v1 + (u * big_stepsize * gradlogtarget(x1) + diagonal_term)*exp_integral_2_b
  x2 <- x1 - v2*(1+4*splitting_const)*big_stepsize
  
  v3 <- exp_2_b*v2 + (u * big_stepsize * gradlogtarget(x2) + diagonal_term)*exp_integral_2_b
  x3 <- x2 + v3*(1+2*splitting_const)*big_stepsize
  
  gradposition <- gradlogtarget(x3)
  
  v4 <- exp_1_b*v3 + (u * big_stepsize * gradposition + diagonal_term)*exp_integral_1_b
  
  return(list(position = x3,
              momentum = v4 - sigma*(H_big_increment - 6*K_big_increment),
              gradposition = gradposition))
}

# Kernel for coupled ULD chain
coupled_kernel <- function(chains) {
  small_step_chain <- chains$small_step_chain
  big_step_chain <- chains$big_step_chain
  
  W_big_increment <- zero_vec
  H_big_increment <- zero_vec
  K_big_increment <- zero_vec
  
  small_time <- 0
  
  for (i in 1:no_of_small_steps) {
    W_small_increment <- t(fast_rmvnorm(1, mean = zero_vec, covariance = W_small_var))
    H_small_increment <- t(fast_rmvnorm(1, mean = zero_vec, covariance = H_small_var))
    K_small_increment <- t(fast_rmvnorm(1, mean = zero_vec, covariance = K_small_var))

    
    small_step_chain <- ULD_small_step(small_step_chain, W_small_increment, H_small_increment, K_small_increment)
    
    
    K_big_increment <- K_big_increment + small_stepsize_squared *
                        ((third*W_small_increment) + 0.5*H_small_increment - K_small_increment) +
                        small_stepsize * (small_stepsize_halfed*W_big_increment + small_time * 
                        (W_big_increment + 0.5*W_small_increment + H_small_increment))
    
    H_big_increment <- H_big_increment + small_stepsize *
                        (W_big_increment + 0.5*W_small_increment + H_small_increment)
    
    W_big_increment <- W_big_increment + W_small_increment
    
    small_time <- small_time + small_stepsize
  }
  
  # Compute the time areas for the Brownian path over the course increment
  H_big_increment <- H_big_increment/big_stepsize - 0.5*W_big_increment
  
  K_big_increment <- third*W_big_increment + 0.5*H_big_increment -
                      K_big_increment/big_stepsize_squared
  
  big_step_chain <- ULD_big_step(big_step_chain, W_big_increment, H_big_increment, K_big_increment)
  
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