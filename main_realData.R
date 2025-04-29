###############################################################################
# Topic: Bayesian Econometrics
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Replicate Conley, Hansen and Rossi (2012)'s paper ``Plausibly Exogenous'' for Acemoglu, Johnson and Robinson (2001) data
###############################################################################
# Organize the working environment
###############################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI", "MASS", "MCMCpack", "matrixsampling", "progress", "foreign")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
remove(install.lib, lib, load.lib)

# Set the seed
set.seed(20250415)

# Import the functions
source("functions.R")  # Import the functions

###############################################################################
# Load the data
###############################################################################

# Load the data on data/acemoglu.dta file
raw_data <- read.dta("data/maketable5.dta")  # Load the data

# Select the relevant columns and rename them
ds <- raw_data[c("shortnam","logpgp95", "avexpr", "logem4")] 
ds$W <- 1
colnames(ds) <- c("ID","Y", "X", "Z", "W")  # Rename the columns
ds <- na.omit(ds)  # Remove missing values

# Remove the following countries from the dataset as Table A1 in the paper
# CHN, KOR, THA, SUR, GBR, FRA
ds <- ds[!(ds$ID %in% c("CHN", "KOR", "THA", "SUR", "GBR", "FRA")), ]  # Remove the countries from the dataset

# Run the MCMC sampler
mu_X = mean(ds$X)  # Mean of X
mu_Y = mean(ds$Y)  # Mean of Y
results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = c(0,mu_X), V_tau = diag(c(1, mu_X**2/4)),
                                     mu_theta = c(0,0,mu_Y), V_theta = diag(c(1, 0.0001, mu_Y**2/4)),
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

###########################################################################
# Diagnostics
###########################################################################

# Plot the trace plots for tau and theta
gg <- ggplot(tau_samples, aes(x = 1:nrow(tau_samples), y = tau_samples[, 1])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(pi))
gg_template_save(gg, "figures/acemoglu/trace_pi.png")

gg <- ggplot(tau_samples, aes(x = 1:nrow(tau_samples), y = tau_samples[, 2])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(phi))
gg_template_save(gg, "figures/acemoglu/trace_phi.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 1])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(beta))
gg_template_save(gg, "figures/acemoglu/trace_beta.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 2])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(gamma))
gg_template_save(gg, "figures/acemoglu/trace_gamma.png")

gg <- ggplot(theta_samples, aes(x = 1:nrow(theta_samples), y = theta_samples[, 3])) + geom_line(size=1.5) + labs(x="Iteration", y=expression(delta))
gg_template_save(gg, "figures/acemoglu/trace_delta.png")

# Histograms for tau and theta
acf <- acf(tau_samples[, 1], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) + 
  geom_bar(stat="identity") + labs(x="Lag", y=expression(pi)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/acemoglu/acf_pi.png")

acf <- acf(tau_samples[, 2], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(phi)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/acemoglu/acf_phi.png")

acf <- acf(theta_samples[, 1], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(beta)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/acemoglu/acf_beta.png")

acf <- acf(theta_samples[, 2], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(gamma)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/acemoglu/acf_gamma.png")

acf <- acf(theta_samples[, 3], plot = FALSE)
gg <- ggplot(data.frame(lag = acf$lag, acf = acf$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(delta)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg_template_save(gg, "figures/acemoglu/acf_delta.png")

# Histograms for tau and theta
gg <- ggplot(data.frame(val = tau_samples[, 1]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(pi), y = "Density") + 
gg_template_save(gg, "figures/acemoglu/hist_pi.png")

gg <- ggplot(data.frame(val = tau_samples[, 2]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(phi), y = "Density") + 
gg_template_save(gg, "figures/acemoglu/hist_phi.png")

gg <- ggplot(data.frame(val = theta_samples[, 1]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(beta), y = "Density") + 
gg_template_save(gg, "figures/acemoglu/hist_beta.png")

gg <- ggplot(data.frame(val = theta_samples[, 2]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(gamma), y = "Density") + 
gg_template_save(gg, "figures/acemoglu/hist_gamma.png")

gg <- ggplot(data.frame(val = theta_samples[, 3]), aes(x = val)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  labs(x = expression(delta), y = "Density") + 
gg_template_save(gg, "figures/acemoglu/hist_delta.png")

# Create tex table for effective sample size
replication_results <- data.frame(
  "Parameter" = c("pi", "phi", "beta", "gamma", "delta"),
  "True Value" = c(-0.61, "-", 0.94, 0, "-"),
  "Mean" = c(mean(tau_samples[, 1]), mean(tau_samples[, 2]), mean(theta_samples[, 1]), mean(theta_samples[, 2]), mean(theta_samples[, 3])),
  "Standard Deviation" = c(sd(tau_samples[, 1]), sd(tau_samples[, 2]), sd(theta_samples[, 1]), sd(theta_samples[, 2]), sd(theta_samples[, 3])),
  "0.025 Quantile" = c(quantile(tau_samples[, 1], 0.025), quantile(tau_samples[, 2], 0.025), quantile(theta_samples[, 1], 0.025), quantile(theta_samples[, 2], 0.025), quantile(theta_samples[, 3], 0.025)),
  "0.975 Quantile" = c(quantile(tau_samples[, 1], 0.975), quantile(tau_samples[, 2], 0.975), quantile(theta_samples[, 1], 0.975), quantile(theta_samples[, 2], 0.975), quantile(theta_samples[, 3], 0.975)),
  "Effective Sample Size" = c(effectiveSize(tau_samples[, 1]), effectiveSize(tau_samples[, 2]), effectiveSize(theta_samples[, 1]), effectiveSize(theta_samples[, 2]), effectiveSize(theta_samples[, 3]))
)
colnames(replication_results) <- c("Parameter", "True Value", "Mean", "Standard Deviation", "0.025", "0.975", "Effective Sample Size")  # Set column

stargazer(replication_results, type = "latex", title = "Replication Results", summary = FALSE, rownames = FALSE, digits.extra = 0, digits = 2, out = "tables/acemoglu/replication_results.tex", label = "tab:replication_results")


##############################################################################
# Iterate over the variance of gamma    
##############################################################################

# Initialize a vector with the true values of gamma
gamma_variances <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100)  # Variances of gamma

# Initialize a matrix to store the means of theta
# Empty matrix with 0 rows and 5 columns
theta_stats <- matrix(0, nrow = 0, ncol = 5)  # Initialize an empty matrix to store the results
colnames(theta_stats) <- c("Gamma Variance","beta_Mean", "beta_Variance", "beta_Lower_CI", "beta_Upper_CI")  # Set column names for the matrix

# Set the progress bar
current_theta_stats <- matrix(0, nrow = length(gamma_variances), ncol = 5)  # Initialize a matrix to store the results for the current gamma value
for (j in 1:length(gamma_variances)) {    
    # Run the MCMC sampler for each variance of gamma
    results <- MCMC_plausibly_exogenous(ds, n_iter = 10000, burn_in = 2000,
                                     mu_tau = c(0,mu_X), V_tau = diag(c(1, mu_X**2/4)),
                                     mu_theta = c(0,0,mu_Y), V_theta = diag(c(1, gamma_variances[j], mu_Y**2/4)),
                                        nu0 = 3, S0 = diag(2))  # Run the MCMC sampler

    # Extract the samples
    mu_theta_post <- results$mu_theta_post 
    V_theta_post <- results$V_theta_post


    # Calculate the posterior mean, variance, and 95% credible interval for beta
    current_theta_stats[j, 1] <- gamma_variances[j]  # Variance of gamma
    current_theta_stats[j, 2] <- mu_theta_post[1] # Posterior mean of beta
    current_theta_stats[j, 3] <- V_theta_post[1, 1]  # Variance of beta
    current_theta_stats[j, 4] <- mu_theta_post[1] - 1.96 * sqrt(V_theta_post[1, 1])  # Lower CI for beta
    current_theta_stats[j, 5] <- mu_theta_post[1] + 1.96 * sqrt(V_theta_post[1, 1])  # Upper CI for beta
}
theta_stats <- rbind(theta_stats, current_theta_stats)  # Combine the results

#################################################################################
# Plot the results
#################################################################################

# Convert the matrix to a data frame for plotting
theta_stats <- as.data.frame(theta_stats)  # Convert the matrix to a data frame
colnames(theta_stats) <- c("Prior_Variance", "Beta_Mean", "Beta_Variance", "Beta_Lower_CI", "Beta_Upper_CI")  # Set column names for the data frame

# Create the plot
gg <- ggplot(theta_stats, aes(x = log10(Prior_Variance), y = Beta_Mean)) +
  geom_line(size = 1.5) +  # Line for the posterior mean
  geom_ribbon(aes(ymin = Beta_Lower_CI, ymax = Beta_Upper_CI), alpha = 0.2, color = NA) +  # Confidence intervals
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
  filename = "figures/acemoglu/plot_prior_gamma.png",
  plot = gg,
  width = 12,
  height = 8
)
