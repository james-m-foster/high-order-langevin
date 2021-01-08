library(unbiasedmcmc)

set.seed(1)

# Load dataset
data(germancredit)
X <- scale(X)
X[,1] <- rep(1, nrow(X))

n <- nrow(X)
p <- ncol(X)
design_matrix <- unname(X)
tdesign_matrix <- t(design_matrix)
response <- Y
new_response <- 2*response - 1  # map to {-1,1}

nsamples <- nrow(design_matrix)
dimension <- ncol(design_matrix)

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

# Log sigmoid function
stable_log_sigmoid <- function(x){
  output <- vector(mode = "logical", length = length(x))
  mask <- x > 0
  nmask <- !mask
  output[mask] <- -log(1+exp(-x[mask]))
  output[nmask] <- x[nmask] - log1p(exp(x[nmask]))
  return(output)
}

# Sigmoid function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Log density of the posterior (up to a constant)
logtarget <- function(beta){
  xbeta <- design_matrix %*% beta
  loglikelihood <- sum(stable_log_sigmoid(new_response*xbeta))
  return(loglikelihood - sum(beta ** 2)/(2*sigma2))
}

# Gradient of log density of the posterior
gradlogtarget <- function(beta){
  xbeta <- design_matrix %*% beta
  tdesign_matrix %*% (sigmoid(-new_response * xbeta) * new_response) - beta / sigma2
}
