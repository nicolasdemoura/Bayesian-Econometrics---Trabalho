###############################################################################
# Topic: Bayesian Econometrics
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Replicate Conley, Hansen and Rossi (2012)'s paper ``Plausibly Exogenous'' for simulated data
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
# First Study: Main Simulation
###############################################################################

# Set parameters for the simulation
N <- 150  # Number of observations

# Coefficients for the structural equation and first stage equation
beta_true <- 2
gamma_true <- 0
delta_true <- -1
pi_true <- 2
phi_true <- 1
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix for the errors

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(10, 10)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(10, 0.0001, 10)),
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
gg <- ggplot(tau_samples, aes(x = 1:nrow(tau_samples), y = tau_samples[, 1])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(pi))
gg_template_save(gg, "figures/simulation/trace_pi.png")

gg <- ggplot(tau_samples, aes(x = 1:nrow(tau_samples), y = tau_samples[, 2])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(phi))
gg_template_save(gg, "figures/simulation/trace_phi.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 1])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(beta))
gg_template_save(gg, "figures/simulation/trace_beta.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 2])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(gamma))
gg_template_save(gg, "figures/simulation/trace_gamma.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 3])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(delta))
gg_template_save(gg, "figures/simulation/trace_delta.png")

# ACF for tau and theta
acf <- acf(tau_samples[, 1], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) + 
  geom_bar(stat="identity") + labs(x="Lag", y=expression(pi)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/simulation/acf_pi.png")

acf <- acf(tau_samples[, 2], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(phi)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/simulation/acf_phi.png")

acf <- acf(theta_samples[, 1], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(beta)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/simulation/acf_beta.png")

acf <- acf(theta_samples[, 2], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(gamma)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/simulation/acf_gamma.png")

acf <- acf(theta_samples[, 3], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(delta)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/simulation/acf_delta.png")

# Histograms for tau and theta
gg <- ggplot(data.frame(val = tau_samples[, 1]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(pi), y = "Density") + 
  geom_vline(xintercept = pi_true, color = "red", linetype = "dashed", size = 1)
gg_template_save(gg, "figures/simulation/hist_pi.png")

gg <- ggplot(data.frame(val = tau_samples[, 2]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(phi), y = "Density") + 
  geom_vline(xintercept = phi_true, color = "red", linetype = "dashed", size = 1)
gg_template_save(gg, "figures/simulation/hist_phi.png")

gg <- ggplot(data.frame(val = theta_samples[, 1]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(beta), y = "Density") + 
  geom_vline(xintercept = beta_true, color = "red", linetype = "dashed", size = 1)
gg_template_save(gg, "figures/simulation/hist_beta.png")

gg <- ggplot(data.frame(val = theta_samples[, 2]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(gamma), y = "Density") + 
  geom_vline(xintercept = gamma_true, color = "red", linetype = "dashed", size = 1)
gg_template_save(gg, "figures/simulation/hist_gamma.png")

gg <- ggplot(data.frame(val = theta_samples[, 3]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(delta), y = "Density") + 
  geom_vline(xintercept = delta_true, color = "red", linetype = "dashed", size = 1)
gg_template_save(gg, "figures/simulation/hist_delta.png")

# Create tex table for effective sample size
simulation_results <- data.frame(
  "Parameter" = c("pi", "phi", "beta", "gamma", "delta"),
  "True Value" = c(pi_true, phi_true, beta_true, gamma_true, delta_true),
  "Mean" = c(mean(tau_samples[, 1]), mean(tau_samples[, 2]), mean(theta_samples[, 1]), mean(theta_samples[, 2]), mean(theta_samples[, 3])),
  "Standard Deviation" = c(sd(tau_samples[, 1]), sd(tau_samples[, 2]), sd(theta_samples[, 1]), sd(theta_samples[, 2]), sd(theta_samples[, 3])),
  "0.025 Quantile" = c(quantile(tau_samples[, 1], 0.025), quantile(tau_samples[, 2], 0.025), quantile(theta_samples[, 1], 0.025), quantile(theta_samples[, 2], 0.025), quantile(theta_samples[, 3], 0.025)),
  "0.975 Quantile" = c(quantile(tau_samples[, 1], 0.975), quantile(tau_samples[, 2], 0.975), quantile(theta_samples[, 1], 0.975), quantile(theta_samples[, 2], 0.975), quantile(theta_samples[, 3], 0.975)),
  "Effective Sample Size" = c(effectiveSize(tau_samples[, 1]), effectiveSize(tau_samples[, 2]), effectiveSize(theta_samples[, 1]), effectiveSize(theta_samples[, 2]), effectiveSize(theta_samples[, 3]))
)
colnames(simulation_results) <- c("Parameter", "True Value", "Mean", "Standard Deviation", "0.025", "0.975", "Effective Sample Size")  # Set column

stargazer(simulation_results, type = "latex", title = "Simulation Results", summary = FALSE, rownames = FALSE, digits.extra = 0, digits = 2, out = "tables/simulation/simulation_results.tex", label = "tab:simulation_results")


##############################################################################
# Second Study: Four different models 
##############################################################################

# Set parameters for the simulation
N <- 150  # Number of observations

########### Baseline model 

# Coefficients for the structural equation and first stage equation
beta_true <- 2
gamma_true <- 0
delta_true <- -1
pi_true <- 2
phi_true <- 1
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Covariance matrix for the errors

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(10, 10)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(10, 0.0001, 10)),
                                     nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

# Extract the samples
theta_samples <- results$theta_samples  # Extract theta samples

model <- c(
  mean(theta_samples[, 1]),  # Posterior mean of beta
  sd(theta_samples[, 1]),  # Posterior standard deviation of beta
  quantile(theta_samples[, 1], 0.025),  # Lower CI for beta
  quantile(theta_samples[, 1], 0.975)   # Upper CI for beta
)
models <- as.data.frame(model)  # Combine the results

########### Weak Instrument model -> Add noise to Z

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data
ds$Z <- ds$Z + rnorm(N, 0, 10)  # Add noise to Z

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(10, 10)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(10, 0.0001, 10)),
                                     nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

# Extract the samples
theta_samples <- results$theta_samples  # Extract theta samples

model <- c(
  mean(theta_samples[, 1]),  # Posterior mean of beta
  sd(theta_samples[, 1]),  # Posterior standard deviation of beta
  quantile(theta_samples[, 1], 0.025),  # Lower CI for beta
  quantile(theta_samples[, 1], 0.975)   # Upper CI for beta
)
models <- cbind(models, model)  # Combine the results

########### Uncertain Instrument model -> Increase prior variance of gamma

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(10, 10)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(10, 1, 10)),
                                     nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

# Extract the samples
theta_samples <- results$theta_samples  # Extract theta samples

model <- c(
  mean(theta_samples[, 1]),  # Posterior mean of beta
  sd(theta_samples[, 1]),  # Posterior standard deviation of beta
  quantile(theta_samples[, 1], 0.025),  # Lower CI for beta
  quantile(theta_samples[, 1], 0.975)   # Upper CI for beta
)
models <- cbind(models, model)  # Combine the results

########### Invalid Instrument model -> gamma_true != 0

gamma_true <- 10

# Simulate data
ds <- simulate_data(N, beta_true, gamma_true, delta_true, pi_true, phi_true, Sigma)  # Simulate data

# Run the MCMC sampler
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = rep(0, 2), V_tau = diag(c(10, 10)),
                                     mu_theta = rep(0, 3), V_theta = diag(c(10, 0.0001, 10)),
                                     nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

# Extract the samples
theta_samples <- results$theta_samples  # Extract theta samples

model <- c(
  mean(theta_samples[, 1]),  # Posterior mean of beta
  sd(theta_samples[, 1]),  # Posterior standard deviation of beta
  quantile(theta_samples[, 1], 0.025),  # Lower CI for beta
  quantile(theta_samples[, 1], 0.975)   # Upper CI for beta
)

models <- cbind(models, model)  # Combine the results
models <- t(models)

colnames(models) <- c("Mean", "Standard Deviation", "Lower CI", "Upper CI")  # Set row names for the data frame
rownames(models) <- c("Baseline", "Weak Instrument", "Uncertain Instrument", "Invalid Instrument")  # Set column names for the data frame
models

# Create tex table for models
stargazer(models, type = "latex", title = "Simulation Results", summary = FALSE, rownames = TRUE, digits.extra = 0, digits = 2, out = "tables/simulation/models.tex", label = "tab:models")

##############################################################################
# Third Study: Iterate over the variance of gamma    
##############################################################################

# Initialize a vector with the true values of gamma
gamma_true_values <- c(0, 1, 5, 10, 20)  # True values of gamma
gamma_variances <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)

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
    title = "",
    x = "Log10(Prior Variance of Gamma)",
    y = "Posterior Mean (and CI) of Beta"
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
  filename = "figures/simulation/plot_prior_gamma.png",
  plot = gg,
  width = 12,
  height = 8
)
