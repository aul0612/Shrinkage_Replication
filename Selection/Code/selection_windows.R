rm(list = ls())
#setwd("Replication")
#%%
# Cross-validation to select optimal lambda
library(glmnet)
library(dplyr)
library(cv)
library(MASS)
library(doParallel)
library(foreach)
library(progress)
library(pbapply)
library(stargazer)
library(lme4)

SEED <- 123
set.seed(SEED)

# Function to calculate BIC
calculate_bic <- function(model, X, y) {
  n <- nrow(X)  # Number of observations
  y_pred <- predict(model, X)  # Predictions for all lambdas
  residuals <- sweep(y_pred, 1, y, FUN = "-")  # Calculate residuals
  rss <- colSums(residuals^2)  # Residual sum of squares for each lambda
  k <- model$df  # Number of predictors for each lambda
  tLL <- model$nulldev - deviance(model)  # Total log likelihood
  bic <- log(n) * k - tLL  # BIC formula
  return(bic)
}

# Define function for Lasso + OLS procedure
lasso_ols <- function(data, lambda, k) {
  Y <- data[, 1]
  X <- data[, -1]
  lasso_model <- glmnet(X, Y, lambda = lambda)
  coeffs <- coef(lasso_model)
  nonzero_coeffs_bool <- coeffs[-1] != 0
  n_nonzero_coeffs <- sum(nonzero_coeffs_bool)
  
  true_coeffs_bool <- c(rep(TRUE, 5), rep(FALSE, k - 5))
  check_true_model <- all(nonzero_coeffs_bool == true_coeffs_bool)
  
  X_relevant <- X[, nonzero_coeffs_bool]
  df <- as.data.frame(cbind(Y, X_relevant))
  ols_model <- lm(Y ~ ., data = df)
  
  coef_orig_length <- numeric(k + 1)
  coef_orig_length[1] <- coef(ols_model)[1]
  j <- 2
  for (i in 2:(k + 1)) {
    if (nonzero_coeffs_bool[i - 1]) {
      coef_orig_length[i] <- coef(ols_model)[j]
      j <- j + 1
    } else {
      coef_orig_length[i] <- 0
    }
  }
  return(list(
    ols_model = ols_model,
    check_true_model = check_true_model,
    n_nonzero_coeffs = n_nonzero_coeffs,
    coeffs_ols = coef_orig_length
  ))
}

# Set parameters
n <- 500
k <- 20
k_rel <- 5
k_irr_corr <- 5
k_irr_ucorr <- k - k_irr_corr - k_rel
B <- 1000
sigma <- 0.2
irr_corr_param <- 10
means_X <- round(runif(k, -10, 10), digits = 2)
lambda_cv <- seq(0.001, 0.5, length.out = 1000)
lambda_bic <- seq(0.001, 0.5, length.out = 1000)
lambda_values <- lambda_bic
rho <- 0.3

# Create valid correlation matrix
Sigma <- matrix(0, k, k)
diag(Sigma) <- 4
block_size <- k_irr_corr
for (i in 1:block_size) {
  for (j in 1:block_size) {
    Sigma[i, j + block_size] <- rho
    Sigma[j + block_size, i] <- rho
  }
}
if (min(eigen(Sigma)$values) <= 0) stop("Correlation matrix is not positive definite")

seeds <- sample(0:10000, B)
beta <- c(rnorm(5, mean = 10, sd = 1), rep(0, k - 5))
X_train <- mvrnorm(n, means_X, Sigma)

# Register parallel backend for Windows
cluster <- makeCluster(detectCores() - 2)
clusterExport(cluster, varlist = c("seeds", "n", "k", "k_rel", "k_irr_corr", 
                              "k_irr_ucorr", "irr_corr_param", "means_X", "X_train",
                              "Sigma", "beta", "lambda_cv", "lambda_bic", 
                              "lambda_values", "calculate_bic", "lasso_ols", 
                              "B", "sigma", "SEED", "cv.glmnet", "cv"))
registerDoParallel(cluster)

# Use pblapply for parallel processing with progress bar
results <- pblapply(1:B, function(b) {
  set.seed(seeds[b])
  eps_train <- rnorm(n, 0, sigma)
  Y_train <- X_train %*% beta + eps_train
  benchmark_ols <- lm(Y_train ~ X_train)
  coeffs_bench_all <- coef(benchmark_ols)
  benchmark_ols_true <- lm(Y_train ~ X_train[, 1:k_rel])
  coeffs_ols_true <- as.vector(c(coef(benchmark_ols_true), rep(0, k - k_rel)))
  
  model_obj <- cv.glmnet(X_train, Y_train, parallel = TRUE, lambda = lambda_cv)
  coeffs <- coef(model_obj, s = "lambda.min")[1:(k + 1)]
  check_true_model_cv <- (all(coeffs[2:6] != 0) && all(coeffs[7:(k + 1)] == 0))
  n_nonzero_coeffs_cv <- sum(coeffs[2:(k + 1)] != 0)
  lambda_min_cv <- model_obj$lambda.min
  
  glmnet_fit_bic <- glmnet(X_train, Y_train, alpha = 1, lambda = lambda_bic)
  bic_values <- calculate_bic(glmnet_fit_bic, X_train, Y_train)
  lambda_min_bic <- glmnet_fit_bic$lambda[which.min(bic_values)]
  coeffs_bic <- matrix(coef(glmnet_fit_bic, s = lambda_min_bic))
  check_true_model_bic <- (all(coeffs_bic[2:6] != 0) && all(coeffs_bic[7:(k + 1)] == 0))
  n_nonzero_coeffs_bic <- sum(coeffs_bic[2:(k + 1)] != 0)
  
  cv_criterion_all_ols <- numeric(length(lambda_values))
  data <- cbind(Y_train, as.data.frame(X_train))
  coeffs_ols <- matrix(NA, nrow = length(lambda_values), ncol = k + 1)
  
  for (l in seq_along(lambda_values)) {
    results_lasso_ols <- lasso_ols(data, lambda_values[l], k)
    cv_criterion_all_ols[l] <- suppressMessages(cv(results_lasso_ols$ols_model, seed = SEED, quietly = TRUE)$`CV crit`)
    coeffs_ols[l, ] <- results_lasso_ols$coeffs_ols
  }
  
  lambda_min_idx <- which.min(cv_criterion_all_ols)
  coeffs_all_ols <- coeffs_ols[lambda_min_idx, ]
  
  return_vector <- c(
    lambda_min_cv,
    lambda_min_bic,
    check_true_model_cv,
    check_true_model_bic,
    coeffs,
    coeffs_all_ols,
    coeffs_bic,
    coeffs_bench_all,
    coeffs_ols_true
  )
  return(return_vector)
}, cl = cluster)

# Convert results to matrix and save to file
results <- do.call(rbind, results)
results_summary <- data.frame(
  MSE = c(mean(results[, 1]), mean(results[, 2])),
  Prob_True_Model = c(mean(results[, 3]), mean(results[, 4]))
)
write.csv(results_summary, file = "Selection/data/selection.csv")

stopCluster(cluster)
stargazer(results_summary, summary = FALSE, rownames = TRUE, title = "Comparison of Parameter Selection Procedures", out = "Selection/tables/selection.tex")
#%%