###############################################################################
# Functions 
###############################################################################


###############################################################################
# Simulate Data
###############################################################################

# Structural Equation:
# Y = beta * X + gamma * Z + delta * W + e1
# X: Endogenous variable
# Y: Dependent variable
# Z: Instrumental variable
# W: Exogenous variable

# First stage equation:
# X = pi * Z + phi * W + e2

simulate_data <- function(N, beta, gamma, delta, pi, phi, Sigma) {
    # N: Number of observations
    # beta: Coefficient for X on Y
    # gamma: Coefficient for Z on Y
    # delta: Coefficient for W on Y
    # pi: Coefficient for Z on X
    # phi: Coefficient for W on X
    # Sigma: Covariance matrix for errors

    e <- mvrnorm(N, mu = c(0, 0), Sigma = Sigma)  # Errors
    e1 <- e[, 1]  # Error term for Y
    e2 <- e[, 2]  # Error term for X

    # Generate W as a column of ones
    W <- matrix(1, nrow = N, ncol = 1)  # W is a column of ones

    # Generate Z as a random normal variable
    Z <- matrix(rnorm(N), nrow = N, ncol = 1)  # Z is a random normal variable

    # Generate X using the first stage equation
    X <- Z %*% pi + W %*% phi + e2  

    # Generate y using the structural equation
    Y <- X %*% beta + Z %*% gamma + W %*% delta + e1

    return(data.frame(Y = Y, X = X, Z = Z, W = W))
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
    Y <- ds$Y
    X <- ds$X
    Z <- ds$Z
    W <- ds$W

    # Initialize matrices to store the samples
    tau_samples <- matrix(0, nrow = n_iter, ncol = 2)  # Coefficients from Z, W to X
    theta_samples <- matrix(0, nrow = n_iter, ncol = 3)  # Coefficients for Z, W and X to Y
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
    # Initialize the progress bar
    pb <- progress_bar$new(total = n_iter, format = "  Progress [:bar] :percent eta: :eta", clear = FALSE)

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
        #Calculate X_hat given tau
        X_hat <- Z * tau_current[1] + W * tau_current[2]  # First stage equation
        XZW <- cbind(X_hat, Z, W)  # Combine X, Z and W
        sigma2_theta <- Sigma_current[1, 1]  # Variance of structural error
        # Update the posterior mean and covariance for theta
        V_theta_post <- solve(V_theta_inv + t(XZW) %*% XZW / sigma2_theta)
        mu_theta_post <- V_theta_post %*% (V_theta_inv %*% mu_theta + t(XZW) %*% Y / sigma2_theta)
        # Draw theta from the posterior distribution
        theta_current <- mvrnorm(1, mu = mu_theta_post, Sigma = V_theta_post)  # Draw theta from the posterior distribution

        ### Update Sigma given tau and theta
        # Calculate the residuals
        e1 <- Y - X * theta_current[1] - Z * theta_current[2] - W * theta_current[3]  # Structural equation residuals
        e2 <- X - Z * tau_current[1] - W * tau_current[2]  # First stage equation residuals
        # Calculate the residual sum of squares
        RSS <- t(cbind(e1, e2)) %*% cbind(e1, e2)  # Residual sum of squares
        
        # Update the posterior parameters for Sigma
        nu <- nu0 + N  # Degrees of freedom
        S <- S0 + RSS  # Scale matrix
        # Draw Sigma from the inverse-Wishart distribution
        Sigma_current <- riwish(nu, S)  # Draw Sigma from the inverse-Wishart distribution

        ### Store the samples
        tau_samples[i, ] <- tau_current  # Store tau samples
        theta_samples[i, ] <- theta_current  # Store theta samples
        Sigma_samples[[i]] <- Sigma_current  # Store Sigma samples
    }

    # Remove burn-in samples
    tau_samples <- tau_samples[-(1:burn_in), ]  # Remove burn-in samples for tau
    theta_samples <- theta_samples[-(1:burn_in), ]  # Remove burn-in samples for theta
    Sigma_samples <- Sigma_samples[-(1:burn_in)]  # Remove burn-in samples for Sigma

    # Column names for the samples
    colnames(tau_samples) <- c("pi", "phi")  # Column names for tau samples
    colnames(theta_samples) <- c("beta", "gamma", "delta")  # Column names for theta samples

    # Get the Sigma_hat from the average of the samples
    Sigma_hat <- apply(simplify2array(Sigma_samples), c(1, 2), mean)  # Posterior mean of Sigma
    sigma2_tau_hat <- Sigma_hat[2, 2]  # Posterior mean of tau
    sigma2_theta_hat <- Sigma_hat[1, 1]  # Posterior mean of theta

    # Get the posterior for tau and theta
    ZW <- cbind(Z, W)  # Combine W and Z
    # Update the posterior mean and covariance for tau
    V_tau_post <- solve(V_tau_inv + t(ZW) %*% ZW / sigma2_tau_hat)  
    mu_tau_post <- V_tau_post %*% (V_tau_inv %*% mu_tau + t(ZW) %*% X / sigma2_tau_hat) 
    # Draw tau from the posterior distribution

    ### Update theta given Sigma
    XZW <- cbind(X, Z, W)  # Combine X, Z and W
    # Update the posterior mean and covariance for theta
    V_theta_post <- solve(V_theta_inv + t(XZW) %*% XZW / sigma2_theta_hat)
    mu_theta_post <- V_theta_post %*% (V_theta_inv %*% mu_theta + t(XZW) %*% Y / sigma2_theta_hat)
    

    return(list(tau_samples = tau_samples, theta_samples = theta_samples, Sigma_samples = Sigma_samples,
                Sigma_hat = Sigma_hat, 
                mu_tau_post = mu_tau_post, V_tau_post = V_tau_post,
                mu_theta_post = mu_theta_post, V_theta_post = V_theta_post))
}

#################################################################################
# Plotting function
#################################################################################

# Plot template function
gg_template_save <- function(p, filename) {
  p +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
          legend.position = "bottom",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 16),
          legend.key.size = unit(1, "cm"),
          legend.background = element_rect(color = "black", size = 0.5)) -> p2
  ggsave(filename, plot = p2, width = 10, height = 10)
}

# Create table for acceptance rates and effective sample size
effectiveSize <- function(x) {
  n <- length(x)
  acf_x <- acf(x, plot = FALSE)$acf[-1]
  var_x <- var(x)
  ess <- n / (1 + 2 * sum(acf_x))
  return(ess)
}