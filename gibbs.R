###############################################################################
# Topic: Learning Gibbs Sampler
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Simulate a Gibbs Sampler for a simple linear model with a normal prior 
###############################################################################
# Organize the working environment
###############################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron")
# Load necessary libraries
library(MASS)
set.seed(123)

# Simulate instruments and errors
n <- 500
Z <- matrix(rnorm(n), n, 1)
true_pi <- 1.5
true_beta <- 2.0

# Simulate correlated errors
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
errors <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
nu <- errors[, 1]
eps <- errors[, 2]

# Generate D and Y
D <- Z %*% true_pi + nu
Y <- D * true_beta + eps

# Initialize
n_iter <- 5000
beta_draws <- numeric(n_iter)
pi_draws <- numeric(n_iter)
Sigma_draws <- array(NA, dim = c(2, 2, n_iter))

# Starting values
beta <- 0
pi <- 0
Sigma <- diag(2)
mu_pi_prior <- 0  # Prior mean for pi
V_pi_prior <- diag(2) * 0.01  # Prior covariance for pi

for (s in 1:n_iter) {
  # Residuals
  nu <- D - Z %*% pi
  eps <- Y - D * beta
  
  # Sample Sigma (2x2 covariance matrix) from inverse-Wishart
  S <- cbind(nu, eps)
  df <- n + 2  # prior df + n
  scale_matrix <- t(S) %*% S
  Sigma <- solve(rWishart(1, df, solve(scale_matrix))[,,1])
  
  # Conditional means and variances
  Sigma11 <- Sigma[1,1]
  Sigma12 <- Sigma[1,2]
  Sigma22 <- Sigma[2,2]
  Sigma21 <- Sigma[2,1]
  
  # Sample pi | rest
  V_pi <- solve(solve(V_pi_prior) + t(Z) %*% Sigma %*% Z)
  mu_pi <- V_pi %*% (solve(V_pi_prior) %*% mu_pi_prior + t(Z) %*% Sigma %*% (Z))
  pi <- rnorm(1, mean = mu_pi, sd = sqrt(V_pi))
  
  # Sample beta | rest
  V_beta <- 1 / sum(D^2) * Sigma22
  mu_beta <- sum(D * Y) / sum(D^2)
  beta <- rnorm(1, mean = mu_beta, sd = sqrt(V_beta))
  
  # Store draws
  pi_draws[s] <- pi
  beta_draws[s] <- beta
  Sigma_draws[,,s] <- Sigma
}

burn_in <- 1000
hist(beta_draws[burn_in:n_iter], breaks = 40,
     main = expression(paste("Posterior of ", beta)),
     xlab = expression(beta))
abline(v = true_beta, col = "red", lwd = 2)

mean(beta_draws[burn_in:n_iter])
quantile(beta_draws[burn_in:n_iter], c(0.025, 0.975))
