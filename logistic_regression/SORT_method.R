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
sixth <- 1/6
two_thirds <- 2/3

exp_1_s <- exp(-0.5*gamma*small_stepsize)
exp_2_s <- exp(-gamma*small_stepsize)
exp_integral_1_s <- (1 - exp_1_s)/(gamma)
exp_integral_2_s <- (1 - exp_2_s)/(gamma)
exp_integral_3_s <- (exp_1_s + 0.5*gamma*small_stepsize - 1)/(gamma ** 2)
exp_integral_4_s <- (exp_1_s + 0.5*gamma*small_stepsize - 1)/(small_stepsize * (gamma ** 2))
exp_integral_5_s <- (exp_2_s + gamma*small_stepsize - 1)/(gamma ** 2)
exp_integral_6_s <- (exp_2_s + gamma*small_stepsize - 1)/(small_stepsize * (gamma ** 2))
exp_integral_7_s <- (1 - exp_2_s)/(small_stepsize*gamma)

exp_1_b <- exp(-0.5*gamma*big_stepsize)
exp_2_b <- exp(-gamma*big_stepsize)
exp_integral_1_b <- (1 - exp_1_b)/(gamma)
exp_integral_2_b <- (1 - exp_2_b)/(gamma)
exp_integral_3_b <- (exp_1_b + 0.5*gamma*big_stepsize - 1)/(gamma ** 2)
exp_integral_4_b <- (exp_1_b + 0.5*gamma*big_stepsize - 1)/(big_stepsize * (gamma ** 2))
exp_integral_5_b <- (exp_2_b + gamma*big_stepsize - 1)/(gamma ** 2)
exp_integral_6_b <- (exp_2_b + gamma*big_stepsize - 1)/(big_stepsize * (gamma ** 2))
exp_integral_7_b <- (1 - exp_2_b)/(big_stepsize*gamma)

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

# Small step of the SORT method
SORT_small_step <- function(chain_state, W_small_increment, H_small_increment, K_small_increment) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum + sigma*(H_small_increment + 6*K_small_increment)
  
  gradposition1 <- chain_state$gradposition
  
  diagonal_term <- sigma*(W_small_increment - 12*K_small_increment)
  
  x1 <- x0 + exp_integral_1_s*v0 + exp_integral_3_s*u*gradposition1 + exp_integral_4_s*diagonal_term
  
  gradposition2 <- gradlogtarget(x1)
  
  x2 <- x0 + exp_integral_2_s*v0 + exp_integral_5_s*u*(third*gradposition1 + two_thirds*gradposition2)  +
                                   exp_integral_6_s * diagonal_term
  
  gradposition3 <- gradlogtarget(x2)
  
  v1 <- exp_2_s*v0 + u*(sixth*exp_2_s*gradposition1 + two_thirds*exp_1_s*gradposition2 +
                        sixth*gradposition3) * small_stepsize +
                        exp_integral_7_s * diagonal_term
  
 
  return(list(position = x2,
              momentum = v1 - sigma*(H_small_increment - 6*K_small_increment),
              gradposition = gradposition3))
}

# Big step of the SORT method
SORT_big_step <- function(chain_state, W_big_increment, H_big_increment, K_big_increment) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum + sigma*(H_big_increment + 6*K_big_increment)
  
  gradposition1 <- chain_state$gradposition
  
  diagonal_term <- sigma*(W_big_increment - 12*K_big_increment)
  
  x1 <- x0 + exp_integral_1_b*v0 + exp_integral_3_b*u*gradposition1 + exp_integral_4_b*diagonal_term
  
  gradposition2 <- gradlogtarget(x1)
  
  x2 <- x0 + exp_integral_2_b*v0 + exp_integral_5_b*u*(third*gradposition1 + two_thirds*gradposition2)  +
                                   exp_integral_6_b * diagonal_term
  
  gradposition3 <- gradlogtarget(x2)
  
  v1 <- exp_2_b*v0 + u*(sixth*exp_2_b*gradposition1 + two_thirds*exp_1_b*gradposition2 +
                        sixth*gradposition3) * big_stepsize +
                        exp_integral_7_b * diagonal_term
  
  
  return(list(position = x2,
              momentum = v1 - sigma*(H_big_increment - 6*K_big_increment),
              gradposition = gradposition3))
  
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

    
    small_step_chain <- SORT_small_step(small_step_chain, W_small_increment, H_small_increment, K_small_increment)
    
    
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
  
  big_step_chain <- SORT_big_step(big_step_chain, W_big_increment, H_big_increment, K_big_increment)
  
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
