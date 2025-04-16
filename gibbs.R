rm(list = ls())  # Clear the workspace

library(MASS)  # For mvrnorm function
library(dplyr)  # For data manipulation
library(MCMCpack)
library(matrixsampling)
# Progress bar
library(progress)  # For progress bar

# Gibbs sampler step for beta
gibbs_sampler_beta <- function(y, X, Z, mu_beta_prior, V_beta_prior, mu_pi_prior, V_pi_prior, s_prior, nu_prior, n_iter = 1000, burn_in = 100) {
  
  # Storage matrices
  n_beta <- length(mu_beta_prior)
  n_pi <- length(mu_pi_prior)
  beta_draws <- matrix(NA, nrow = n_iter, ncol = n_beta)
  pi_draws <- matrix(NA, nrow = n_iter, ncol = n_pi)
  H_draws <- array(NA, dim = c(nrow(X), nrow(X), n_iter))
  H_beta_draws <- array(NA, dim = c(nrow(X), nrow(X), n_iter))
  H_pi_draws <- array(NA, dim = c(nrow(Z), nrow(Z), n_iter))

  # Set initial values
  beta_draws[1, ] <- mu_beta_prior 
  pi_draws[1, ] <- mu_pi_prior
  H_draws[, , 1] <- diag(100,nrow(X))  
  H_beta_draws[, , 1] <- diag(100,nrow(X))  
  H_pi_draws[, , 1] <- diag(100,nrow(Z))

  # Precompute inverses
  V_beta_prior_inv <- solve(V_beta_prior)
  V_pi_prior_inv <- solve(V_pi_prior)
  
  # Set up progress bar
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent eta: :eta",
    total = n_iter - 1,
    clear = FALSE,
    width = 60
  )
  iter <- 2  # Initialize iteration counter
  for (iter in 2:n_iter) {
    # Update progress bar
    pb$tick()

    # Select 2x2 submatrix of H for each equation
    H_pi_draws[,,iter-1] <- H_draws[1:2, 1:2, iter-1]  # First stage equation
    H_beta_draws[,,iter-1] <- H_draws[3:4, 3:4, iter-1]  # Structural equation

    # Sample from posterior for Pi conditional on H and Beta
    V_pi_n <- solve(V_pi_prior_inv + t(Z) %*% H_pi_draws[,,iter-1] %*% Z)
    mu_pi_n <- V_pi_n %*% (V_pi_prior_inv %*% mu_pi_prior + t(Z) %*% H_pi_draws[,,iter-1] %*% y)
    pi_sample <- MASS::mvrnorm(1, mu = mu_pi_n, Sigma = V_pi_n)
    
    # Sample from posterior for Beta conditional on H and Pi
    V_beta_n <- solve(V_beta_prior_inv + t(X) %*% H_beta_draws[,,iter-1] %*% X)
    mu_beta_n <- V_beta_n %*% (V_beta_prior_inv %*% mu_beta_prior + t(X) %*% H_beta_draws[,,iter-1] %*% y)
    beta_sample <- MASS::mvrnorm(1, mu = mu_beta_n, Sigma = V_beta_n)

    # Store draw
    pi_draws[iter, ] <- pi_sample
    beta_draws[iter, ] <- beta_sample

    # Sample from posterior for H conditional on Beta and Pi
    eps <- y - X %*% beta_sample
    u <- X[,2] - Z %*% pi_sample 
    eps_tilde <- cbind(eps, u) 

    # Sample from the inverse-Wishart distribution
    s_n <- s_prior + eps_tilde %*% t(eps_tilde)  # Posterior scale parameter
    nu_n <- nu_prior + nrow(X)  # Posterior degrees of freedom
    H_sample <- rinvwishart(1, nu_n, s_n)  # Inverse-Wishart sample

    # Store draw
    H_draws[, , iter] <- H_sample  # Store draw
  }

  # Drop burn-in samples
  beta_tilde_draws_burnin <- beta_tilde_draws[(burn_in + 1):n_iter, , drop = FALSE]
  H_draws_burnin <- H_draws[, , (burn_in + 1):n_iter, drop = FALSE]

  return(list(
    beta_draws = beta_tilde_draws_burnin,
    H_draws = H_draws_burnin,
    all_beta_draws = beta_tilde_draws,
    all_H_draws = H_draws
  ))
}

# Create a function to simulate data for a simple IV model
simulate_data <- function(n, beta, pi, Sigma) {
  # Get two correlated error terms, one for the structural equation and one for the first stage
  e <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  
  # Generate the Instrumental variable (Z) 
  Z <- cbind(1, rnorm(n))  

  # Generate the endogenous variable (X) using the first stage equation
  X <- Z %*% pi + e[, 1] 
  X <- cbind(1, X)  # Add intercept to X

  # Generate the dependent variable (Y) using the structural equation
  y <- X %*% beta + e[, 2]

  return(data.frame(y = y, X = X, Z = Z))
}

n <- 50  # Number of observations
beta <- as.matrix(c(-15,2), transpose = TRUE)  # Coefficient for the structural equation
pi <- as.matrix(c(20, 2.5), transpose = TRUE)  # Coefficient for the first stage equation
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix for the error terms

ds <- simulate_data(n, beta, pi, Sigma)  # Simulate data 
X <- cbind(ds$X.1, ds$X.2)  # Design matrix (including intercept)
Z <- cbind(ds$Z.1, ds$Z.2)  # Instrumental variable matrix (including intercept)
y <- ds$y  # Dependent variable

# Set the parameters for the Gibbs sampler
mu_beta_prior <- c(0,0)
V_beta_prior <- diag(c(50,1))
mu_pi_prior <- c(0,0)
V_pi_prior <- diag(c(50,1))

s_prior <- diag(5, nrow = nrow(X))  # Prior scale matrix for the error terms
nu_prior <- 5  # Prior degrees of freedom for the error terms

# Set the n_iter and burn_in parameters
n_iter <- 10000  # Number of iterations for the Gibbs sampler
burn_in <- 1000  # Number of burn-in iterations

# Run the Gibbs sampler
gibbs_results <- gibbs_sampler_beta(y, X, Z, mu_beta_prior, V_beta_prior, mu_pi_prior, V_pi_prior, s_prior, nu_prior, n_iter, burn_in) 

# Extract the posterior draws
beta_draws <- gibbs_results$beta_draws

# Plot the time series of the posterior draws for beta
par(mfrow = c(2, 2))
plot(beta_draws[, 1], type = "l", main = "Posterior Draws for Beta_0", ylab = "Value", xlab = "Iteration")
plot(beta_draws[, 2], type = "l", main = "Posterior Draws for Beta_1", ylab = "Value", xlab = "Iteration")
plot(beta_draws[, 3], type = "l", main = "Posterior Draws for Beta_2", ylab = "Value", xlab = "Iteration")
plot(beta_draws[, 4], type = "l", main = "Posterior Draws for Beta_3", ylab = "Value", xlab = "Iteration")

# Print the summary statistics of the posterior draws
summary_stats <- apply(beta_draws, 2, function(x) c(mean = mean(x), sd = sd(x), quantiles = quantile(x, probs = c(0.025, 0.5, 0.975))))
print(summary_stats)
