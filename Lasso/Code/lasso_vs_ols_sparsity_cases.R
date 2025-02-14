rm(list = ls())

library(parallel)
library(glmnet) 
library(ggplot2)  
library(MASS)       
library(tidyr)
library(dplyr)

# Function to generate data and perform Monte Carlo simulations
monte_carlo_simulation <- function(n, p, n_sim, true_coefficients, lambdas, correlation, seed = NULL) {
  # Set the seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create covariance matrix for correlated predictors
  cov_matrix <- matrix(correlation, nrow = p, ncol = p)
  diag(cov_matrix) <- 1  # Set diagonal elements to 1 (variance of 1)
  
  # Store results for OLS and Lasso
  results <- data.frame()
  
  simulate_once <- function(seed, lambda) {
    set.seed(seed)
    
    # Generate synthetic data with correlated predictors
    X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_matrix)
    epsilon <- rnorm(n, 0, 2)
    y <- X %*% true_coefficients + epsilon
    
    # OLS regression
    ols_fit <- lm(y ~ X - 1)
    ols_coef <- coef(ols_fit)
    
    # Lasso regression
    lasso_fit <- glmnet(X, y, alpha = 1, lambda = lambda, standardize = TRUE, intercept = FALSE)
    lasso_coef <- as.vector(coef(lasso_fit, s = lambda))[-1]  # Exclude intercept
    
    # Calculate statistics for both methods
    list(
      OLS = list(coef = ols_coef),
      Lasso = list(coef = lasso_coef)
    )
  }
  
  for (lambda in lambdas) {
    # Run simulations in parallel
    n_cores <- detectCores() - 1
    seeds <- 1:n_sim
    cluster <- makeCluster(n_cores)
    clusterExport(cluster, c("simulate_once", "n", "p", "true_coefficients", "lambda", "cov_matrix"), envir = environment())
    clusterEvalQ(cluster, library(glmnet))
    clusterEvalQ(cluster, library(MASS))
    
    simulation_results <- parLapply(cluster, seeds, function(seed) simulate_once(seed, lambda))
    stopCluster(cluster)
    
    # Process results
    for (method in c("OLS", "Lasso")) {
      all_coefs <- do.call(rbind, lapply(simulation_results, function(res) res[[method]]$coef))
      bias <- colMeans(all_coefs) - true_coefficients
      variance <- apply(all_coefs, 2, var)
      mse <- rowMeans((all_coefs - matrix(true_coefficients, nrow = n_sim, ncol = p, byrow = TRUE))^2)
      
      results <- rbind(results, data.frame(
        Lambda = lambda,
        Method = method,
        Bias = mean(bias^2),
        Variance = mean(variance),
        MSE = mean(mse)
      ))
    }
  }
  
  return(results)
}

# Parameters for the simulation
n <- 100             # Number of observations
p <- 15              # Number of predictors
n_sim <- 100         # Number of Monte Carlo simulations
lambdas <- seq(0.00, 0.5, by = 0.02)  # Range of lambda values
correlation <- 0     # Correlation coefficient for predictors
seed <- 1234           # Seed for reproducibility

# Varying sparsity levels
sparsity_levels <- list(
  `Dense (15 non-zero)` = c(runif(15, -5, 5)),
  `Moderate (10 non-zero)` = c(runif(10, -5, 5), rep(0, 5)),
  `Sparse (5 non-zero)` = c(runif(5, -5, 5), rep(0, 10))
)

# Run simulation for each sparsity level
all_results <- data.frame()
for (sparsity_name in names(sparsity_levels)) {
  true_coefficients <- sparsity_levels[[sparsity_name]]
  results <- monte_carlo_simulation(n, p, n_sim, true_coefficients, lambdas, correlation, seed)
  results$Sparsity <- sparsity_name
  all_results <- rbind(all_results, results)
}

# Reshape and plot
plot_data <- all_results %>%
  pivot_longer(cols = c(Bias, Variance, MSE), names_to = "Metric", values_to = "Value")
 
custom_titles <- c('Dense (15 non-zero)'="(a)", 'Moderate (10 non-zero)'= "(b)",'Sparse (5 non-zero)'="(c)")

p <- ggplot(plot_data, aes(x = Lambda, y = Value, color = Metric, linetype = Method)) +
  geom_line() +
  facet_wrap(~ Sparsity, labeller = labeller(Sparsity = custom_titles)) +  # Facet by sparsity level
  theme_minimal() +
  geom_line(linewidth = 0.7) +
  labs(y = "Metric Value", x = expression(lambda)) +
  scale_color_manual(values = c("Bias" = "#e81717", "MSE" = "#000000", "Variance" = "#2417e8")) +
  theme(legend.position = "bottom")

ggsave("Lasso/Output/lasso_sparsity_cases.png",
       width = 25, height = 16, units = "cm",
       plot = p,
       device = cairo_pdf)