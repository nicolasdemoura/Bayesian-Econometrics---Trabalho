    ###############################################################################
    # Functions 
    ###############################################################################

    ###############################################################################
    # Prior, Likelihood, and Posterior
    ###############################################################################

    # Prior
    log_prior <- function(delta, beta, gamma1, gamma2, Sigma) {
        

        return(lp)  # Prior up to constant
    }

    # Log-likelihood
    log_likelihood <- function(y, tau, beta, lambda, sigma2) {
        ll <- sum(dnorm(y, mean = nelson_siegel(tau, beta, lambda), sd = sqrt(sigma2), log = TRUE))
        return(ll)  # Likelihood up to constant
    }

    # Log-posterior
    log_posterior <- function(y, tau, beta, lambda, sigma2) {
        return(log_prior(beta, lambda, sigma2) + log_likelihood(y, tau, beta, lambda, sigma2))
    }


###############################################################################
# Simulate Data
###############################################################################

simulate_data <- function(N, beta, gamma1, gamma2, pi, sigma2) {
    # Generate the error terms
    e <- rnorm(N, 0, sigma2)
    v <- rnorm(N, 0, sigma2)

    # Generate the instrumental variable
    Z <- matrix(rnorm(N * length(pi)), nrow = N, ncol = length(pi))

    # Generate the exogenous variable
    W <- matrix(rnorm(N * length(gamma1)), nrow = N, ncol = length(gamma1))

    # Generate the endogenous variable
    X <- matrix(rnorm(N * length(beta)), nrow = N, ncol = length(beta))
    X <- Z %*% t(pi) + v + e

    # Generate the outcome variable
    Y <- X %*% t(beta) + W %*% t(gamma1) + Z %*% t(gamma2) + e

    # Return the data as a data frame
    data <- data.frame(Y = Y, X = X, Z = Z, W = W)
    # Set w as a matrix
    data$W <- as.matrix(data$W)
    data$Z <- as.matrix(data$Z)
    data$X <- as.vector(data$X)
    data$Y <- as.vector(data$Y)
    colnames(data) <- c("Y", paste0("X", 1:length(beta)), paste0("Z", 1:length(pi)), paste0("W", 1:length(gamma1)))
    return(data)
}

###############################################################################
# Plot the data
###############################################################################

plot_and_save <- function(data1, data2, filename, x_col, y_col, color_col, plot_type = "curve") {
    library(ggplot2)
    
    # Combine the datasets
    combined_data <- rbind(data1, data2)
    
    # Create the plot based on the plot_type argument
    if (plot_type == "curve") {
        p <- ggplot(combined_data, aes_string(x = x_col, y = y_col, color = color_col)) +
            geom_line(size = 1.5)
    } else if (plot_type == "scatter") {
        p <- ggplot(combined_data, aes_string(x = x_col, y = y_col, color = color_col)) +
            geom_point(size = 3) +
            geom_smooth(method = "lm", se = FALSE, size = 1.5)
    } else {
        stop("Invalid plot_type. Use 'curve' or 'scatter'.")
    }
    
    # Add common plot elements
    p <- p +
        labs(x = x_col,
             y = y_col,
             color = color_col) +
        theme_bw(base_size = 25) +
        theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
              legend.position = "bottom",
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 16),
              legend.key.size = unit(1, "cm"),
              legend.background = element_rect(color = "black", size = 0.5)) 
                  
    # Save the plot
    ggsave(filename, plot = p, width = 10, height = 10)
}
