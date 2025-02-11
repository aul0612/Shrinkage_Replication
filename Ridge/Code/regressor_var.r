rm(list = ls())
#setwd("Replication")
library(MASS)
library(parallel)
library(Matrix)
library(progress)
library(caret)
library(ggplot2)
library(dplyr)

# set seed
set.seed(374)


# set parameters
n <- 50       # sample size
k <- 20         # number of regressors
B  <- 1000   # number of simulations
sigma  <- 2 # std of errors

# set some values for the true parameter vector
beta <- rep(1, k)
means_X <- rep(0, k)
vars_X <- seq(.1, 5, length.out=k)
# set up grid of lambda values
lambda_values <- c(seq(0,20, length.out = 100), seq(20, 12000, length.out = 401)[-1])
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

# set up containers
bias2_all <- matrix(NA, nrow = L, ncol = k)
variance_all <- matrix(NA, nrow = L, ncol = k)
coefficients_all <- array(NA, c(L, B, k))
prediction_error_all <- array(NA, c(L, B, n))
mspe_all <- numeric(L)
mean_coefficient <- matrix(NA, L, B)

cov_mat <- matrix(0, k, k)
diag(cov_mat) <- vars_X
X_train <- mvrnorm(n, means_X, cov_mat)
X_test <- mvrnorm(n, means_X, cov_mat)


pb <- progress_bar$new(total = B)
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
        prediction_error_all[i, b, ] <- X_test %*% coefficients_all[i, b, ] - Y_test
    }
}

for (i in 1:L){
    bias2_all[i, ] <- (colMeans(coefficients_all[i,,])-beta)^2
    variance_all[i, ] <- apply(coefficients_all[i,,], 2, var)
    mspe_all[i] <- mean((prediction_error_all[i,,])^2)
    mean_coefficient[i,] <- colMeans(coefficients_all[i,,])
}


# plot the on mean coefficients for on lambda
matplot(log(lambda_values), mean_coefficient, type = "l",
 xlab = expression("log("~lambda~ ")"), ylab = "Ridge Coefficient", main = "Regularization paths for different variances") 

vars_X 
cols <- c("red", "green", "blue", "purple", "orange")
upper <- 1:5
var_df <- data.frame(Var = vars_X, str = NA, col = NA, lab = NA, variable = 1:k)
for (i in 1:length(vars_X)){
  for (j in upper){
    if (vars_X[i]>j-1 && vars_X[i]<=j){
      var_df$lab[i] <- j
      var_df$str[i] <- paste0("σₓ² ∈ (",j-1,",", j, "]")
      var_df$col[i] <- cols[j]
    }
  }
}
var_df

pdf("Ridge/Output/Ridge_Variance.pdf")
matplot(log(lambda_values), mean_coefficient, type = "l",
col = var_df$col,
xlab = expression("log("~lambda~ ")"), ylab = "Ridge Coefficient",
main = "Regularization paths for different variances") 
legend("topright", legend = unique(var_df$str), fill = unique(var_df$col))
dev.off()