###############################################################################
# Topic: Bayesian Econometrics
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Replicate Conley, Hansen and Rossi (2012)'s paper ``Plausibly Exogenous'' 
###############################################################################
# Organize the working environment
###############################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI", "MASS", "MCMCpack", "matrixsampling", "progress")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
remove(install.lib, lib, load.lib)

# Set the seed
set.seed(20250415)

###############################################################################
# Simulate Data
###############################################################################

# Structural Equation:
# y = beta * X + gamma * Z + delta * W + e1

# First stage equation:
# X = pi * Z + phi * W + e2

simulate_data <- function(N, beta, gamma, delta, pi, phi, Sigma) {
    # N: Number of observations
    # beta: Coefficient for X on y
    # gamma: Coefficient for Z on y
    # delta: Coefficient for W on y
    # pi: Coefficient for Z on X
    # phi: Coefficient for W on X
    # Sigma: Covariance matrix for errors

    e <- mvrnorm(N, mu = c(0, 0), Sigma = Sigma)  # Errors
    e1 <- e[, 1]  # Error term for y
    e2 <- e[, 2]  # Error term for X

    # Generate W as a column of ones
    W <- matrix(1, nrow = N, ncol = 1)  # W is a column of ones

    # Generate Z as a random normal variable
    Z <- matrix(rnorm(N), nrow = N, ncol = 1)  # Z is a random normal variable

    # Generate X using the first stage equation
    X <- Z %*% pi + W %*% phi + e2  

    # Generate y using the structural equation
    y <- X %*% beta + Z %*% gamma + W %*% delta + e1

    return(data.frame(y = y, X = X, Z = Z, W = W))
}

###############################################################################
# Gibbs Sampler for Bayesian IV Estimation
###############################################################################

MCMC_plausibly_exogenous <- function(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau =  rep(0, 2), V_tau = diag(c(1, 100)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(1, 1, 100)),
                                     nu0 = 3, S0 = diag(2)) {

    # Get the dimensions of the data
    N <- nrow(ds)  # Number of observations

    # Get the data
    y <- ds$y
    X <- ds$X
    Z <- ds$Z
    W <- ds$W

    # Initialize matrices to store the samples
    tau_samples <- matrix(0, nrow = n_iter, ncol = 2)  # Coefficients from Z, W to X
    theta_samples <- matrix(0, nrow = n_iter, ncol = 3)  # Coefficients for Z, W and X to y
    Sigma_samples <- vector("list", n_iter)

    # Initialize the parameters
    tau_current <- mu_tau  # Initial values for tau
    theta_current <- mu_theta  # Initial values for theta
    Sigma_current <- S0  # Initial value for Sigma

    # Store the first sample
    tau_samples[1, ] <- tau_current
    theta_samples[1, ] <- theta_current
    Sigma_samples[[1]] <- Sigma_current

    # Invert the prior covariance matrices
    V_tau_inv <- solve(V_tau)  # Inverse of the prior covariance matrix for tau
    V_theta_inv <- solve(V_theta)  # Inverse of the prior covariance matrix for theta

    # Set the progress bar
    pb <- progress_bar$new(total = n_iter, format = "  Progress [:bar] :percent eta: :eta", clear = FALSE)
    pb$tick(0)  # Initialize the progress bar

    # Gibbs sampling loop
    for(i in 2:n_iter) {
        pb$tick()  # Update the progress bar
        
        ### Update tau given Sigma
        ZW <- cbind(Z, W)  # Combine W and Z
        sigma2_tau <- Sigma_current[2, 2]  # Variance of first stage error
        # Update the posterior mean and covariance for tau
        V_tau_post <- solve(V_tau_inv + t(ZW) %*% ZW / sigma2_tau)  
        mu_tau_post <- V_tau_post %*% (V_tau_inv %*% mu_tau + t(ZW) %*% X / sigma2_tau) 
        # Draw tau from the posterior distribution
        tau_current <- mvrnorm(1, mu = mu_tau_post, Sigma = V_tau_post)  # Draw tau from the posterior distribution

        ### Update theta given Sigma
        XZW <- cbind(X, Z, W)  # Combine X, Z and W
        sigma2_theta <- Sigma_current[1, 1]  # Variance of structural error
        # Update the posterior mean and covariance for theta
        V_theta_post <- solve(V_theta_inv + t(XZW) %*% XZW / sigma2_theta)
        mu_theta_post <- V_theta_post %*% (V_theta_inv %*% mu_theta + t(XZW) %*% y / sigma2_theta)
        # Draw theta from the posterior distribution
        theta_current <- mvrnorm(1, mu = mu_theta_post, Sigma = V_theta_post)  # Draw theta from the posterior distribution

        ### Update Sigma given tau and theta
        # Calculate the residuals
        e1 <- y - X * theta_current[1] - Z * theta_current[2] - W * theta_current[3]  # Structural equation residuals
        e2 <- X - Z * tau_current[1] - W * tau_current[2]  # First stage equation residuals
        # Calculate the residual sum of squares
        RSS <- t(cbind(e1, e2)) %*% cbind(e1, e2)  # Residual sum of squares
        # Update the posterior parameters for Sigma
        nu0 <- nu0 + N  # Degrees of freedom
        S0 <- S0 + RSS  # Scale matrix
        # Draw Sigma from the inverse-Wishart distribution
        Sigma_current <- riwish(nu0, S0)  # Draw Sigma from the inverse-Wishart distribution

        ### Store the samples
        tau_samples[i, ] <- tau_current  # Store tau samples
        theta_samples[i, ] <- theta_current  # Store theta samples
        Sigma_samples[[i]] <- Sigma_current  # Store Sigma samples
    }

    # Remove burn-in samples
    tau_samples <- tau_samples[-(1:burn_in), ]  # Remove burn-in samples for tau
    theta_samples <- theta_samples[-(1:burn_in), ]  # Remove burn-in samples for theta
    Sigma_samples <- Sigma_samples[-(1:burn_in)]  # Remove burn-in samples for Sigma

    return(list(tau_samples = tau_samples, theta_samples = theta_samples, Sigma_samples = Sigma_samples))
}


###############################################################################
# Main Simulation
###############################################################################

# Set parameters for the simulation
N <- 150  # Number of observations

# Coefficients for the structural equation and first stage equation
beta_true <- 2
gamma_true <- 0
delta_true <- 0
pi_true <- 2
phi_true <- 0

# Covariance matrix for the errors
Sigma <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)  # Covariance matrix for the errors

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(1, 1000)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(1, 0.0001, 1000)),
                                     nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

# Extract the samples
tau_samples <- results$tau_samples  # Extract tau samples
theta_samples <- results$theta_samples  # Extract theta samples
Sigma_samples <- results$Sigma_samples  # Extract Sigma samples

# Column names for the samples
colnames(tau_samples) <- c("pi", "phi")  # Column names for tau samples
colnames(theta_samples) <- c("beta", "gamma", "delta")  # Column names for theta samples

# Summary statistics 
summary(tau_samples)  # Summary statistics for tau samples
summary(theta_samples)  # Summary statistics for theta samples

###########################################################################
# Diagnostics
###########################################################################

# Plot the trace plots for tau and theta
plot(tau_samples[, 1], type = "l", main = "Trace plot for tau[1]", ylab = "tau[1]")  # Trace plot for tau[1]
plot(tau_samples[, 2], type = "l", main = "Trace plot for tau[2]", ylab = "tau[2]")  # Trace plot for tau[2]
plot(theta_samples[, 1], type = "l", main = "Trace plot for theta[1]", ylab = "theta[1]")  # Trace plot for theta[1]
plot(theta_samples[, 2], type = "l", main = "Trace plot for theta[2]", ylab = "theta[2]")  # Trace plot for theta[2]
plot(theta_samples[, 3], type = "l", main = "Trace plot for theta[3]", ylab = "theta[3]")  # Trace plot for theta[3]

# Histograms for tau and theta
hist(tau_samples[, 1], main = "Histogram of tau[1]", xlab = "tau[1]", breaks = 50)  # Histogram of tau[1]
hist(tau_samples[, 2], main = "Histogram of tau[2]", xlab = "tau[2]", breaks = 50)  # Histogram of tau[2]
hist(theta_samples[, 1], main = "Histogram of theta[1]", xlab = "theta[1]", breaks = 50)  # Histogram of theta[1]
hist(theta_samples[, 2], main = "Histogram of theta[2]", xlab = "theta[2]", breaks = 50)  # Histogram of theta[2]
hist(theta_samples[, 3], main = "Histogram of theta[3]", xlab = "theta[3]", breaks = 50)  # Histogram of theta[3]
