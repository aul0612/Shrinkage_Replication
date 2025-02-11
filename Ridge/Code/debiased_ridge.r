

#%%
rm(list = ls())
#setwd("Replication")
library(MASS)
library(parallel)
library(Matrix)
library(progress)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
#%%

# set seed
set.seed(374)


# set parameters
n <- 50      # sample size
k <- 30        # number of regressors
B  <- 1000   # number of simulations
sigma  <- 2 # std of errors

# set some values for the true parameter vector
beta <- rnorm(k, mean=1, sd=0.1)
means_X <- round(runif(k, -10, 10), digits=2)
# set up grid of lambda values
lambda_values <- seq(0, 30, length.out = 100)
L <- length(lambda_values)

# to avoid centering and scaling:
# first consider X values that are mean zero and have the same variance
# and are not correlated
# Also consider model where intercept is 0

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
debiased_ridge_func <- function(X, Y, lambda){
    ridge_coeffs <- ridge_func(X,Y,lambda)
    debiased_ridge_coeffs <- ridge_coeffs + lambda * solve(t(X)%*%X + lambda*diag(ncol(X)))%*%ridge_coeffs
    return(debiased_ridge_coeffs)
}

# set up containers
bias2_all <- matrix(NA, nrow = L, ncol = k)
variance_all <- matrix(NA, nrow = L, ncol = k)
coefficients_all <- array(NA, c(L, B, k))
prediction_error_all <- array(NA, c(L, B, n))
mspe_all <- numeric(L)

# containers for debiased estimator
bias2_all_d <- matrix(NA, nrow = L, ncol = k)
variance_all_d <- matrix(NA, nrow = L, ncol = k)
coefficients_all_d <- array(NA, c(L, B, k))
prediction_error_all_d <- array(NA, c(L, B, n))
mspe_all_d <- numeric(L)

X_train <- mvrnorm(n, means_X, diag(k))
X_test <- mvrnorm(n, means_X, diag(k))
pb <- progress_bar$new(total=B)
for (b in 1:B){
    pb$tick()
    # Simulate Data
    #X_train <- mvrnorm(n, means_X, diag(k))
    #X_test <- mvrnorm(n, means_X, diag(k))
    eps_train <- rnorm(n, 0, sigma)
    eps_test <- rnorm(n, 0, sigma)
    Y_train <- X_train %*%beta + eps_train
    Y_test <- X_test %*% beta + eps_test
    for (i in 1:L){
        coefficients_all[i, b, ] <- ridge_func(X_train, Y_train, lambda_values[i])
        #prediction_error_all[i, b, ] <- X_test %*% coefficients_all[i, b, ] - Y_test
        coefficients_all_d[i, b, ] <- debiased_ridge_func(X_train, Y_train, lambda_values[i])
        #prediction_error_all_d[i, b, ] <- X_test %*% coefficients_all_d[i, b, ] - Y_test
    }
}

for (i in 1:L){
    bias2_all[i, ] <- (colMeans(coefficients_all[i,,])-beta)^2
    variance_all[i, ] <- apply(coefficients_all[i,,], 2, var)
    #mspe_all[i] <- mean((prediction_error_all[i,,])^2)
    bias2_all_d[i, ] <- (colMeans(coefficients_all_d[i,,])-beta)^2
    variance_all_d[i, ] <- apply(coefficients_all_d[i,,], 2, var)
    #mspe_all_d[i] <- mean((prediction_error_all_d[i,,])^2)
}

# Calculate total MSE
bias2_total <- rowSums(bias2_all)
variance_total <- rowSums(variance_all)


mse_total <- bias2_total + variance_total

# Calculate total MSE for debiased
bias2_total_d <- rowSums(bias2_all_d)
variance_total_d <- rowSums(variance_all_d)


mse_total_d <- bias2_total_d + variance_total_d


library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming you have your data in these vectors, first create a dataframe
plot_data_small_sample <- data.frame(
  lambda = lambda_values,
  MSE = mse_total,
  MSE_d = mse_total_d,
  Bias2 = bias2_total,
  Bias2_d = bias2_total_d,
  Variance = variance_total,
  Variance_d = variance_total_d
) %>%
  pivot_longer(cols = -lambda, names_to = "metric", values_to = "value") %>%
  mutate(
    Estimator = ifelse(grepl("_d$", metric), "Debiased Ridge", "Normal Ridge"),
    Metric = gsub("_d$", "", metric)
  )

write.csv(plot_data_small_sample, "Ridge/Output/deb_small_sample.csv")
#%%


#%%   Generate larger sample comparison
# set seed
set.seed(374)


# set parameters
n <- 200      # sample size
k <- 30        # number of regressors
B  <- 1000   # number of simulations
sigma  <- 2 # std of errors

# set some values for the true parameter vector
beta <- rnorm(k, mean=1, sd=0.1)
means_X <- round(runif(k, -10, 10), digits=2)
# set up grid of lambda values
lambda_values <- seq(0, 30, length.out = 100)
L <- length(lambda_values)

# to avoid centering and scaling:
# first consider X values that are mean zero and have the same variance
# and are not correlated
# Also consider model where intercept is 0

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
debiased_ridge_func <- function(X, Y, lambda){
    ridge_coeffs <- ridge_func(X,Y,lambda)
    debiased_ridge_coeffs <- ridge_coeffs + lambda * solve(t(X)%*%X + lambda*diag(ncol(X)))%*%ridge_coeffs
    return(debiased_ridge_coeffs)
}

# set up containers
bias2_all <- matrix(NA, nrow = L, ncol = k)
variance_all <- matrix(NA, nrow = L, ncol = k)
coefficients_all <- array(NA, c(L, B, k))
prediction_error_all <- array(NA, c(L, B, n))
mspe_all <- numeric(L)

# containers for debiased estimator
bias2_all_d <- matrix(NA, nrow = L, ncol = k)
variance_all_d <- matrix(NA, nrow = L, ncol = k)
coefficients_all_d <- array(NA, c(L, B, k))
prediction_error_all_d <- array(NA, c(L, B, n))
mspe_all_d <- numeric(L)

X_train <- mvrnorm(n, means_X, diag(k))
X_test <- mvrnorm(n, means_X, diag(k))
pb <- progress_bar$new(total=B)
for (b in 1:B){
    pb$tick()
    # Simulate Data
    #X_train <- mvrnorm(n, means_X, diag(k))
    #X_test <- mvrnorm(n, means_X, diag(k))
    eps_train <- rnorm(n, 0, sigma)
    #eps_test <- rnorm(n, 0, sigma)
    Y_train <- X_train %*%beta + eps_train
    #Y_test <- X_test %*% beta + eps_test
    for (i in 1:L){
        coefficients_all[i, b, ] <- ridge_func(X_train, Y_train, lambda_values[i])
        #prediction_error_all[i, b, ] <- X_test %*% coefficients_all[i, b, ] - Y_test
        coefficients_all_d[i, b, ] <- debiased_ridge_func(X_train, Y_train, lambda_values[i])
        #prediction_error_all_d[i, b, ] <- X_test %*% coefficients_all_d[i, b, ] - Y_test
    }
}

for (i in 1:L){
    bias2_all[i, ] <- (colMeans(coefficients_all[i,,])-beta)^2
    variance_all[i, ] <- apply(coefficients_all[i,,], 2, var)
    #mspe_all[i] <- mean((prediction_error_all[i,,])^2)
    bias2_all_d[i, ] <- (colMeans(coefficients_all_d[i,,])-beta)^2
    variance_all_d[i, ] <- apply(coefficients_all_d[i,,], 2, var)
    #mspe_all_d[i] <- mean((prediction_error_all_d[i,,])^2)
}

# Calculate total MSE
bias2_total <- rowSums(bias2_all)
variance_total <- rowSums(variance_all)


mse_total <- bias2_total + variance_total

# Calculate total MSE for debiased
bias2_total_d <- rowSums(bias2_all_d)
variance_total_d <- rowSums(variance_all_d)


mse_total_d <- bias2_total_d + variance_total_d

plot_data_large_sample <- data.frame(
  lambda = lambda_values,
  MSE = mse_total,
  MSE_d = mse_total_d,
  Bias2 = bias2_total,
  Bias2_d = bias2_total_d,
  Variance = variance_total,
  Variance_d = variance_total_d
) %>%
  pivot_longer(cols = -lambda, names_to = "metric", values_to = "value") %>%
  mutate(
    Estimator = ifelse(grepl("_d$", metric), "Debiased Ridge", "Normal Ridge"),
    Metric = gsub("_d$", "", metric)
  )

write.csv(plot_data_large_sample, "Ridge/Output/deb_large_sample.csv")
#%%


# Combine the two datasets with a sample size indicator
plot_data <- bind_rows(
  mutate(plot_data_large_sample, Sample_Size = "Large Sample"),
  mutate(plot_data_small_sample, Sample_Size = "Small Sample")
)

# Use factor() to specify the order
plot_data$Sample_Size <- factor(plot_data$Sample_Size, 
                                levels = c("Small Sample", "Large Sample"))

p <- ggplot(plot_data, aes(x = lambda, y = value, color = Metric, linetype = Estimator)) +
  geom_line() +
  scale_color_manual(values = c("MSE" = "black", "Bias2" = "red", "Variance" = "blue")) +
  scale_linetype_manual(values = c("Normal Ridge" = "solid", "Debiased Ridge" = "dashed")) +
  facet_wrap(~ Sample_Size, scales = "free_y") +
  labs(
    x = expression(lambda),
    y = "MSE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(title = "Metric"),
    linetype = guide_legend(title = "Estimator")
  )
p

ggsave("Ridge/Output/debiased.pdf", width = 16, height = 7, units = "cm")

