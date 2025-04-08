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
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
remove(install.lib, lib, load.lib)

# Set the seed
set.seed(20250406)

# Import the functions
source("functions.R")
source("rivGibbspe_suff_sbeta.R")

###############################################################################
# Simulate Data
###############################################################################

# Set the parameters for the simulation
N <- 1000 # Number of observations
beta <- 1 # Coefficients for the endogenous variable
gamma1 <- 2 # Coefficients for the exogenous variable
gamma2 <- 0 # Coefficients for the instrument
pi <- 0.5 # Coefficients for the instrument
sigma2 <- 1 # Variance of the error terms

# Apply the function to simulate the data
ds <- simulate_data(N, beta, gamma1, gamma2, pi, sigma2)
colnames(ds) <- c("y", "x", "z", "w")

# Regressions

ivreg(y ~ x + w | z + w, data = ds) # Instrumental variable regression
