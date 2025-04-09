library(MASS)  # For mvrnorm function
library(dplyr)  # For data manipulation
library(MCMCpack)
library(matrixsampling)
# Progress bar
library(progress)  # For progress bar

# Gibbs sampler step for beta
gibbs_sampler_beta <- function(y, X, Z, mu_beta_tilde_prior, V_beta_tilde_prior, s_prior, nu_prior, n_iter = 1000, burn_in = 100) {
  
  # Storage matrices
  n_beta <- length(mu_beta_tilde_prior)
  beta_tilde_draws <- matrix(NA, nrow = n_iter, ncol = n_beta)
  H_draws <- array(NA, dim = c(nrow(X)+nrow(Z), nrow(X)+nrow(Z), n_iter))

  # Set initial values
  beta_tilde_draws[1, ] <- mu_beta_tilde_prior  # Initialize with prior mean
  H_draws[, , 1] <- diag(100,nrow(X)+nrow(Z))  # Initialize with identity matrix

  # Compute X_tilde
  # Add ncol(Z) columns of 0s to X
  X_zero <- cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(Z)))  # Add columns of zeros to X
  # Add ncol(X) columns of 0s to Z
  Z_zero <- cbind(matrix(0, nrow = nrow(Z), ncol = ncol(X)),Z)  # Add columns of zeros to Z

  X_tilde <- rbind(X_zero, Z_zero)
  y_tilde <- c(y, X[,2])

  # Precompute inverses
  V_beta_tilde_prior_inv <- solve(V_beta_tilde_prior)
  
  # Set up progress bar
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent eta: :eta",
    total = n_iter - 1,
    clear = FALSE,
    width = 60
  )

  for (iter in 2:n_iter) {
    # Update progress bar
    pb$tick()
    
    # Sample from posterior for Beta conditional on H
    V_beta_tilde_n <- solve(V_beta_tilde_prior_inv + t(X_tilde) %*% H_draws[,,iter-1] %*% X_tilde)
    mu_beta_tilde_n <- V_beta_tilde_n %*% (V_beta_tilde_prior_inv %*% mu_beta_tilde_prior + t(X_tilde) %*% H_draws[,,iter-1] %*% y_tilde)
    beta_tilde_sample <- MASS::mvrnorm(1, mu = mu_beta_tilde_n, Sigma = V_beta_tilde_n)

    # Store draw
    beta_tilde_draws[iter, ] <- beta_tilde_sample
  
    # Sample from posterior for H conditional on Beta
    eps_tilde <- y_tilde - X_tilde %*% beta_tilde_sample  # Residuals
    s_n <- s_prior + eps_tilde %*% t(eps_tilde)  # Posterior scale parameter
    nu_n <- nu_prior + nrow(X_tilde)  # Posterior degrees of freedom
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

n <- 150  # Number of observations
beta <- as.matrix(c(30,2), transpose = TRUE)  # Coefficient for the structural equation
pi <- as.matrix(c(50, 0.25), transpose = TRUE)  # Coefficient for the first stage equation
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix for the error terms

ds <- simulate_data(n, beta, pi, Sigma)  # Simulate data 
X <- cbind(ds$X.1, ds$X.2)  # Design matrix (including intercept)
Z <- cbind(ds$Z.1, ds$Z.2)  # Instrumental variable matrix (including intercept)
y <- ds$y  # Dependent variable

# Set the parameters for the Gibbs sampler
mu_beta_tilde_prior <- c(40,1,40,1)
V_beta_tilde_prior <- diag(5, nrow = 4)  # Prior covariance matrix for beta
s_prior <- diag(5, nrow = nrow(X)+nrow(Z))  # Prior scale matrix for the error terms
nu_prior <- 5  # Prior degrees of freedom for the error terms

# Set the n_iter and burn_in parameters
n_iter <- 10000  # Number of iterations for the Gibbs sampler
burn_in <- 1000  # Number of burn-in iterations

# Run the Gibbs sampler
gibbs_results <- gibbs_sampler_beta(y, X, Z, mu_beta_tilde_prior, V_beta_tilde_prior, s_prior, nu_prior, n_iter, burn_in) 

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
