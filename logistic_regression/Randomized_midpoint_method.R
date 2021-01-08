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

# Small step of the randomized midpoint method
Midpoint_small_step <- function(chain_state, small_midpoint, small_W1, small_W2, small_W3) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  
  x1 <- x0 + ((1 - exp(-gamma*small_midpoint))/gamma)*v0 + 
              0.5*u*(small_midpoint + ((exp(-gamma*small_midpoint) - 1)/gamma))*gradlogtarget(x0) +
              sigma*small_W1
  
  gradlogtarget_x <- gradlogtarget(x1)
  
  v1 <- exp_s*(v0 + u*small_stepsize*exp(gamma*small_midpoint)*gradlogtarget_x) + sigma*small_W2 
  
  x2 <- x0 + exp_s_gamma*v0 + 0.5*u*small_stepsize*(1-exp(gamma*(small_midpoint - small_stepsize))) *
                                  gradlogtarget_x + sigma*small_W3
  
  return(list(position = x2,
              momentum = v1))
}

# Big step of the randomized midpoint method
Midpoint_big_step <- function(chain_state, big_midpoint, big_W1, big_W2, big_W3) {
  x0 <- chain_state$position
  v0 <- chain_state$momentum
  
  
  x1 <- x0 + ((1 - exp(-gamma*big_midpoint))/gamma)*v0 + 
              0.5*u*(big_midpoint + ((exp(-gamma*big_midpoint) - 1)/gamma))*gradlogtarget(x0) +
              sigma*big_W1
  
  gradlogtarget_x <- gradlogtarget(x1)
  
  v1 <- exp_b*(v0 + u*big_stepsize*exp(gamma*big_midpoint)*gradlogtarget_x) + sigma*big_W2 
  
  x2 <- x0 + exp_b_gamma*v0 + 0.5*u*big_stepsize*(1-exp(gamma*(big_midpoint - big_stepsize))) * 
                                  gradlogtarget_x + sigma*big_W3
  
  return(list(position = x2,
              momentum = v1))
}

# Generate stochatic integrals over the interval [t, t + stepsize]
generate_integrals  <- function(stepsize) {
  smallstepsize <- stepsize
  bigstepsize <- 2*stepsize
  
  exps <- exp(-gamma*smallstepsize)
  expb <- exp(-gamma*bigstepsize)
  
  first_integral <- t(fast_rmvnorm(1, mean = zero_vec, 
                                      covariance = diag((1 - expb)/(2*gamma), 
                                      dimension, dimension)))
  second_integral  <- t(fast_rmvnorm(1, mean = ((1-exps) / (gamma*(1 + exps))) * first_integral,
                                        covariance = diag(((4*exps - expb + 2*gamma*smallstepsize - 3)/
                                                             (2*(gamma^3))) - 
                                                            (((1 - exps)^3) / (2*(gamma^3)*(1+exps))),
                                        dimension, dimension)))
  
  return (list(exp_increment = first_integral, exp_area = second_integral))
}

# Combine stochatic integrals over the intervals [t, t + stepsize1] and [t + stepsize1, t + stepsize2]
combine_integrals <- function(first_integrals, second_integrals, stepsize2) {
  
  integral_1 <- exp(-gamma*stepsize2)*first_integrals$exp_increment + second_integrals$exp_increment
  
  integral_2 <- first_integrals$exp_area + second_integrals$exp_area + 
                   ((1 - exp(-gamma*stepsize2))/gamma)*first_integrals$exp_increment
  
  return (list(exp_increment = integral_1, exp_area = integral_2))
}

# Kernel for coupled ULD chain
coupled_kernel <- function(chains) {
  small_step_chain <- chains$small_step_chain
  big_step_chain <- chains$big_step_chain
  
  midpoint1 <- runif(1)
  midpoint2 <- 1 - midpoint1
  midpoint3 <- runif(1)
  midpoint4 <- 1 - midpoint3
  
  first_integrals <- generate_integrals(midpoint1*small_stepsize)
  second_integrals <- generate_integrals(midpoint2*small_stepsize)
  
  third_integrals <- generate_integrals(midpoint3*small_stepsize)
  fourth_integrals <- generate_integrals(midpoint4*small_stepsize)
  
  first_half_integrals <- combine_integrals(first_integrals, second_integrals,
                                            midpoint2*small_stepsize)
  second_half_integrals <- combine_integrals(third_integrals, fourth_integrals,
                                             midpoint4*small_stepsize)
  
  big_integrals <- combine_integrals(first_half_integrals, second_half_integrals, small_stepsize)
  
  small_step_chain <- Midpoint_small_step(small_step_chain,
                                          midpoint1*small_stepsize,
                                          first_integrals$exp_area,
                                          first_half_integrals$exp_increment,
                                          first_half_integrals$exp_area)
  small_step_chain <- Midpoint_small_step(small_step_chain, 
                                          midpoint3*small_stepsize, 
                                          third_integrals$exp_area, 
                                          second_half_integrals$exp_increment, 
                                          second_half_integrals$exp_area)
  
  # Sample randomized midpoint over the "big" interval
  if (sample(x = 0:1, size = 1) == 0)
  {
    bigmidpoint <- 0.5*midpoint1 
    
    midpoint_integral <- first_integrals$exp_area
  }
  else
  {
    bigmidpoint <- 0.5 + 0.5*midpoint3 
    
    midpoint_integral <- combine_integrals(first_half_integrals, third_integrals, midpoint3*small_stepsize)$exp_area
    
  }

  big_step_chain <- Midpoint_big_step(big_step_chain, bigmidpoint*big_stepsize, midpoint_integral, big_integrals$exp_increment, big_integrals$exp_area)
  
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
