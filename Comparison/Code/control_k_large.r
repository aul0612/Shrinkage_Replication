rm(list = ls())

#%%

library(MASS)
library(parallel)
library(glmnet)
library(mosaic)
library(ggplot2)
library(tidyr)
library(dplyr)
library(progress)

#mu_values <- c(0, 1, 10)
#sigma_values <- c(0.001, 1, 2, 5)
#for (mu in mu_values){
#  for (sigma_beta in sigma_values){
#    print(mu)
#    print(sigma_beta)

    # set seed
    set.seed(374)


    # set parameters
    n <- 50       # sample size
    k <- 20         # number of regressors
    B  <- 1000   # number of simulations
    sigma  <- 2 # std of errors
    #correlation <- 0.6
    #cov_mat <- matrix(correlation, k, k)
    #diag(cov_mat) <- 1

    # set some values for the true parameter vector

    ## beta to show ridge, lasso work fine if true k<n but more parameters involved
    k <- 100
    relevant_k <- 20
    beta <- c(rnorm(relevant_k, 1, 2), rep(0, k-relevant_k))


    means_X <- round(runif(k, -10, 10), digits=2)
    # set up grid of lambda values
    lambda_values <- seq(0, 10, length.out = 100)
    L <- length(lambda_values)


    # define a function that returns the ridge estimator
    ridge_func <- function(X, Y, lambda) {
      A <- t(X) %*% X + (lambda * diag(ncol(X)))
      # handling case of k<n with lambda=0
      if (rankMatrix(A)==k){
        return(solve(t(X) %*% X + (lambda * diag(ncol(X)))) %*% (t(X) %*% Y))
      }
      else{
        return(NA)
      }
    }

    # set up containers
    bias2_all_ridge <- matrix(NA, nrow = L, ncol = k)
    variance_all_ridge <- matrix(NA, nrow = L, ncol = k)
    bias2_all_lasso <- matrix(NA, nrow = L, ncol = k)
    variance_all_lasso <- matrix(NA, nrow = L, ncol = k)
    coefficients_all_ridge <- array(NA, c(L, B, k))
    coefficients_all_lasso <- array(NA, c(L, B, k))
    prediction_error_all_ridge <- array(NA, c(L, B, n))
    prediction_error_all_lasso <- array(NA, c(L, B, n))
    nonzero_coeffs_lasso <- matrix(NA, L, B)
    mspe_all_ridge <- numeric(L)
    mspe_all_lasso <- numeric(L)
    X_train <- mvrnorm(n, means_X, diag(k))
    X_test <- mvrnorm(n, means_X, diag(k))

    #
    pb <- progress_bar$new(total = B)
    for (b in 1:B){
        pb$tick()
        # Simulate Data
        eps_train <- rnorm(n, 0, sigma)
        eps_test <- rnorm(n, 0, sigma)
        Y_train <- X_train %*%beta + eps_train
        Y_test <- X_test %*% beta + eps_test
        for (i in 1:L){
            coefficients_all_ridge[i, b, ] <- ridge_func(X_train, Y_train, lambda_values[i])
            coefficients_all_lasso[i, b, ] <- coef(glmnet(X_train, Y_train, alpha = 1, lambda = lambda_values[i], intercept = FALSE))[-1]
            prediction_error_all_ridge[i, b, ] <- X_test %*% coefficients_all_ridge[i, b, ] - Y_test
            prediction_error_all_lasso[i, b, ] <- X_test %*% coefficients_all_lasso[i, b, ] - Y_test
            nonzero_coeffs_lasso[i, b] <- sum(coefficients_all_lasso[i, b, ]!=0)
        }
    } 

    for (i in 1:L){
        bias2_all_ridge[i, ] <- (colMeans(coefficients_all_ridge[i,,])-beta)^2
        variance_all_ridge[i, ] <- apply(coefficients_all_ridge[i,,], 2, var)
        bias2_all_lasso[i, ] <- (colMeans(coefficients_all_lasso[i,,])-beta)^2
        variance_all_lasso[i, ] <- apply(coefficients_all_lasso[i,,], 2, var)
        mspe_all_ridge[i] <- mean((prediction_error_all_ridge[i,,])^2)
        mspe_all_lasso[i] <- mean((prediction_error_all_lasso[i,,])^2)
    }
    mean_nonzero_coeffs_lasso <- rowMeans(nonzero_coeffs_lasso)
    mspe_ols <- mspe_all_ridge[1]

    bias2_total_ridge <- rowSums(bias2_all_ridge)
    variance_total_ridge <- rowSums(variance_all_ridge)
    mse_total_ridge <- bias2_total_ridge + variance_total_ridge

    bias2_total_lasso <- rowSums(bias2_all_lasso)
    variance_total_lasso <- rowSums(variance_all_lasso)

    mse_total_lasso <- bias2_total_lasso + variance_total_lasso

#%%

data <- data.frame(mean_nonzero_coeffs_lasso,
                  bias2_total_lasso, bias2_total_ridge,
                  variance_total_lasso, variance_total_ridge,
                  mse_total_lasso, mse_total_ridge)

write.csv(data, "Comparison/data/control_k_large.csv")

#%%
