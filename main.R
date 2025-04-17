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

# Import the functions
source("functions.R")  # Import the functions

###############################################################################
# Main Simulation
###############################################################################

# Set parameters for the simulation
N <- 150  # Number of observations

# Coefficients for the structural equation and first stage equation
beta_true <- 2
gamma_true <- 20
delta_true <- 1
pi_true <- 2
phi_true <- 1
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix for the errors

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

##############################################################################
# Iterate over the variance of gamma    
##############################################################################

# Initialize a vector with the true values of gamma
gamma_true_values <- c(0, 1, 5, 10, 20)  # True values of gamma
gamma_variances <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)  # Variances of gamma

# Initialize a matrix to store the means of theta
# Empty matrix with 0 rows and 5 columns
theta_stats <- matrix(0, nrow = 0, ncol = 6)  # Initialize an empty matrix to store the results
rownames(theta_stats) <- paste0("Var_gamma_prior_", gamma_variances)  # Set row names for the matrix
colnames(theta_stats) <- c("True_Gamma","beta_Mean", "beta_Variance", "beta_Lower_CI", "beta_Upper_CI")  # Set column names for the matrix

# Set the progress bar
pb <- progress_bar$new(total = length(gamma_true_values)*length(gamma_variances), format = "  Progress [:bar] :percent eta: :eta", clear = FALSE)

# Loop over the true values of gamma and variances of gamma
for (i in 1:length(gamma_true_values)) {
    # Simulate data for each true value of gamma
    ds <- simulate_data(N, beta_true, gamma_true_values[i], delta_true, pi_true, phi_true, Sigma)  # Simulate data
    
    current_theta_stats <- matrix(0, nrow = length(gamma_variances), ncol = 6)  # Initialize a matrix to store the results for the current gamma value
    for (j in 1:length(gamma_variances)) {    
        pb$tick()  # Update the progress bar

        # Run the MCMC sampler for each variance of gamma
        results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                            mu_tau = rep(0, 2), V_tau = diag(c(1, 100)),
                                            mu_theta = rep(0, 3), V_theta = diag(c(1, gamma_variances[j], 10)),
                                            nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

        # Extract the samples
        mu_theta_post <- results$mu_theta_post 
        V_theta_post <- results$V_theta_post

        # Calculate the posterior mean, variance, and 95% credible interval for beta
        current_theta_stats[j, 1] <- gamma_true_values[i]  # True value of gamma 
        current_theta_stats[j, 2] <- gamma_variances[j]  # Variance of gamma
        current_theta_stats[j, 3] <- mu_theta_post[1] # Posterior mean of beta
        current_theta_stats[j, 4] <- V_theta_post[1, 1]  # Variance of beta
        current_theta_stats[j, 5] <- mu_theta_post[1] - 1.96 * sqrt(V_theta_post[1, 1])  # Lower CI for beta
        current_theta_stats[j, 6] <- mu_theta_post[1] + 1.96 * sqrt(V_theta_post[1, 1])  # Upper CI for beta
    }
    theta_stats <- rbind(theta_stats, current_theta_stats)  # Combine the results
}

#################################################################################
# Plot the results
#################################################################################

# Convert the matrix to a data frame for plotting
theta_stats <- as.data.frame(theta_stats)  # Convert the matrix to a data frame
colnames(theta_stats) <- c("True_Gamma", "Prior_Variance", "Beta_Mean", "Beta_Variance", "Beta_Lower_CI", "Beta_Upper_CI")  # Set column names for the data frame

# Create the plot
gg <- ggplot(theta_stats, aes(x = log10(Prior_Variance), y = Beta_Mean, color = as.factor(True_Gamma), group = True_Gamma)) +
  geom_line(size = 1.5) +  # Line for the posterior mean
  geom_ribbon(aes(ymin = Beta_Lower_CI, ymax = Beta_Upper_CI, fill = as.factor(True_Gamma)), alpha = 0.2, color = NA) +  # Confidence intervals
  scale_color_brewer(palette = "Dark2", name = "True Gamma") +  # Colorblind-friendly palette for lines
  scale_fill_brewer(palette = "Dark2", name = "True Gamma") +  # Colorblind-friendly palette for ribbons
  labs(
    title = "Posterior Mean and CI for Beta by True Gamma",
    x = "Log10(Prior Variance of Gamma)",
    y = "Posterior Mean (and CI)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.margin = unit(c(5, 7, 2, 2), "mm"),
    legend.position = "right",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.background = element_rect(color = "black", size = 0.5)
  )

# Save the plot
ggsave(
  filename = "figures/plot_overlapping_gamma.png",
  plot = gg,
  width = 12,
  height = 8
)
