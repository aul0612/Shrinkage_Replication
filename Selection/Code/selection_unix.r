
rm(list = ls())
#setwd("Replication")
#%%
# Cross validation to select optimal lambda
library(glmnet)
library(dplyr)
library(cv)
library(MASS)
library(doParallel)
library(foreach)
library(doMC)
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
  bic <- log(n)*k - tLL  # BIC formula
  return(bic)
}

# Define function for second procedure (Lasso + OLS)
lasso_ols <- function(data, lambda, k){
    Y <- data[, 1]
    X <- data[,-1]
    lasso_model <- glmnet(X, Y, lambda = lambda)
    coeffs <- coef(lasso_model)
    nonzero_coeffs_bool  <- coeffs[-1] != 0
    n_nonzero_coeffs <- sum(nonzero_coeffs_bool)
    
    # Modify true model check to be dynamic based on k
    true_coeffs_bool <- c(rep(TRUE, 5), rep(FALSE, k-5))
    check_true_model <- all(nonzero_coeffs_bool == true_coeffs_bool)
    
    X_relevant <- X[, nonzero_coeffs_bool]
    df <- as.data.frame(cbind(Y, X_relevant))
    ols_model <- lm(Y ~ ., data = df)
    
    # Dynamic coefficient reconstruction
    coef_orig_length <- numeric(k+1)
    coef_orig_length[1] <- coef(ols_model)[1]
    j <- 2
    for (i in 2:(k+1)){
        if (nonzero_coeffs_bool[i-1] == TRUE){
            coef_orig_length[i] <- coef(ols_model)[j]
            j <- j+1
        }
        else{
            coef_orig_length[i] <- 0
        }
    }
    return(list(ols_model = ols_model,
                check_true_model = check_true_model,
                n_nonzero_coeffs = n_nonzero_coeffs,
                coeffs_ols = coef_orig_length))
}

# Set parameters
n <- 500      # Sample size
k <- 20       # Number of regressors
k_rel <- 5
k_irr_corr <- 5
k_irr_ucorr <- k - k_irr_corr - k_rel
B <- 1000     # Number of simulations
sigma <- .2 # Standard deviation of errors -> 
# .185 works fine for smaller probabilities
# then decreasing to a smaller one to have convergence but also see convergence of both methods while hybrid does not get much better
irr_corr_param <- 10
means_X <- round(runif(k, -10, 10), digits=2)
# Define a fine lambda sequence for BIC
lambda_cv <- seq(0.001, 0.5, length.out=1000)  # Default lambda range for CV
lambda_bic <- seq(0.001, 0.5, length.out = 1000)  # Custom lambda range for BIC
lambda_values <- lambda_bic



rho <- 0.3    # Need to keep this smaller to ensure positive definiteness

# Create valid correlation matrix
Sigma <- matrix(0, k, k)
diag(Sigma) <- 4

# Create block structure
block_size <- k_irr_corr
for(i in 1:block_size) {
    for(j in 1:block_size) {
        # Fill correlation between relevant and corresponding irrelevant variables
        Sigma[i, j+block_size] <- rho
        Sigma[j+block_size, i] <- rho
    }
}

# Verify positive definiteness
eigen_values <- eigen(Sigma)$values
if(min(eigen_values) <= 0) {
    stop("Correlation matrix is not positive definite")
}
# Draw seeds for parallelization loop
seeds <- sample(0:10000,B)

# define relevant lambda range
lambda_values <- seq(0, 1, 0.01)
L <- length(lambda_values)

# Set true parameter vector
beta <- c(rnorm(5, mean=10, sd = 1), rep(0, k-5))

# Simulate X
#X_train_rel <- mvrnorm(n, runif(k_rel, -10, 10), diag(k_rel)) # uncorrelated regressors
#X_train_irr_corr <- X_train_rel + mvrnorm(n, rep(0, k_irr_corr), irr_corr_param*diag(k_irr_corr))
#X_train_irr_ucorr <- mvrnorm(n, runif(k_irr_ucorr, -10, 10), diag(k_irr_ucorr))
#X_train <- cbind(X_train_rel, X_train_irr_corr, X_train_irr_ucorr)


# simulate X with covariance structure Sigma
X_train <- mvrnorm(n, means_X, Sigma)


# Register parallel backend with 4 cores
# Register parallel backend with 4 cores
registerDoParallel(cores = detectCores()-2)

# Use pblapply for parallel processing with progress bar
results <- pblapply(1:B, function(b) {
    set.seed(seeds[b])
    
    # Simulate Data
    eps_train <- rnorm(n, 0, sigma)
    Y_train <- X_train %*% beta + eps_train
    benchmark_ols <- lm(Y_train ~ X_train)
    coeffs_bench_all <- coef(benchmark_ols)
    # benchmark_ols_true
    benchmark_ols_true <- lm(Y_train ~ X_train[, 1:k_rel])
    coeffs_ols_true <- as.vector(c(coef(benchmark_ols_true), rep(0, k-k_rel)))

    # CV for Lasso-CV
    model_obj <- cv.glmnet(X_train, Y_train, parallel = TRUE, lambda = lambda_cv)
    coeffs <- coef(model_obj, s = "lambda.min")[1:(k+1)]

    # Dynamic true model check
    check_true_model_cv <- (all(coeffs[2:6] != 0) && all(coeffs[7:(k+1)] == 0))
    n_nonzero_coeffs_cv <- sum(coeffs[2:(k+1)] != 0)
    lambda_min_cv <- model_obj$lambda.min

    # Fit glmnet model with custom lambda grid for BIC
    glmnet_fit_bic <- glmnet(X_train, Y_train, alpha = 1, lambda = lambda_bic)

    # Calculate BIC for all lambdas in the custom sequence
    bic_values <- calculate_bic(glmnet_fit_bic, X_train, Y_train)
    # Find lambda that minimizes BIC
    lambda_min_bic <- glmnet_fit_bic$lambda[which.min(bic_values)]

    coeffs_bic <- matrix(coef(glmnet_fit_bic, s = lambda_min_bic))
    
    # Dynamic true model check
    check_true_model_bic <- (all(coeffs_bic[2:6] != 0) && all(coeffs_bic[7:(k+1)] == 0))
    n_nonzero_coeffs_bic <- sum(coeffs_bic[2:(k+1)] != 0)

    # CV for Lasso + OLS (Lasso-OLS-CV)
    cv_criterion_all_ols <- numeric(L)
    n_nonzero_coeffs_all_ols <- numeric(L)
    check_true_model_all_ols <- logical(L)

    data <- cbind(Y_train, as.data.frame(X_train))
    coeffs_ols <- matrix(NA, nrow = L, ncol = k+1)

    for (l in 1:L){
        results_lasso_ols <- lasso_ols(data, lambda_values[l], k)
        cv_criterion_all_ols[l] <- suppressMessages(cv(results_lasso_ols$ols_model, seed = SEED, quietly = TRUE)$"CV crit")
        n_nonzero_coeffs_all_ols[l] <- results_lasso_ols$n_nonzero_coeffs
        check_true_model_all_ols[l] <- results_lasso_ols$check_true_model
        coeffs_ols[l, ] <- results_lasso_ols$coeffs_ols
    }

    lambda_min_idx <- which.min(cv_criterion_all_ols)
    lambda_min_ols <- lambda_values[lambda_min_idx]
    n_nonzero_coeffs_ols <- n_nonzero_coeffs_all_ols[lambda_min_idx]
    check_true_model_ols <- check_true_model_all_ols[lambda_min_idx]
    coeffs_all_ols <- coeffs_ols[lambda_min_idx, ]

    
    # Dynamically create return vector based on k
    return_vector <- c(
        lambda_min_cv,
        lambda_min_ols,
        lambda_min_bic,
        check_true_model_cv,
        check_true_model_ols,
        check_true_model_bic,
        n_nonzero_coeffs_cv,
        n_nonzero_coeffs_ols,
        n_nonzero_coeffs_bic
    )
    
    # Append coefficients dynamically
    return_vector <- c(# Alternative correlation formulation## Alternative correlation formulation
        return_vector,
        coeffs,           # coeffs_all_cv
        coeffs_all_ols,   # coeffs_all_ols
        coeffs_bic,       # coeffs_all_bic
        coeffs_bench_all,  # coeffs_all_bench 
        coeffs_ols_true # coeffs_ols_true
    )
    
    return(return_vector)
}, cl = detectCores()-2)

# Convert results to matrix
results <- do.call(rbind, results)

# Dynamic result extraction
# Calculate start indices for each coefficient set
coeff_start_idx <- 10  # After the first 10 summary statistics
cv_coeff_start <- coeff_start_idx
ols_coeff_start <- cv_coeff_start + (k+1)
bic_coeff_start <- ols_coeff_start + (k+1)
bench_coeff_start <- bic_coeff_start + (k+1)
ols_true_coeff_start <- bench_coeff_start + (k+1)

# Extract results dynamically
lambda_min_cv_all <- results[, 1]
lambda_min_ols_all <- results[, 2]
lambda_min_bic_all <- results[, 3]
check_true_model_cv_all <- results[, 4]
check_true_model_ols_all <- results[, 5]
check_true_model_bic_all <- results[, 6]
n_nonzero_coeffs_cv_all <- results[, 7]
n_nonzero_coeffs_ols_all <- results[, 8]
n_nonzero_coeffs_bic_all <- results[, 9]

# Dynamic coefficient extraction
coeffs_all_cv <- results[, cv_coeff_start:(cv_coeff_start+k)]
coeffs_all_ols <- results[, ols_coeff_start:(ols_coeff_start+k)]
coeffs_all_bic <- results[, bic_coeff_start:(bic_coeff_start+k)]
coeffs_all_bench <- results[, bench_coeff_start:(bench_coeff_start+k)]
coeffs_ols_true <- results[, ols_true_coeff_start:(ols_true_coeff_start+k)]



# MSE for OLS-Benchmark
beta_int <- c(0, beta)
bias2_all_bench <- (colMeans(coeffs_all_bench) - beta_int)^2
variance_all_bench <- apply(coeffs_all_bench, MARGIN = 2, var)
mse_all_bench <- bias2_all_bench + variance_all_bench

total_bias2_bench <- sum(bias2_all_bench)
total_variance_bench <- sum(variance_all_bench)
total_mse_bench <- total_bias2_bench + total_variance_bench

# Statistics for Lasso-CV
beta_int <- c(0, beta)
bias2_all_cv <- (colMeans(coeffs_all_cv) - beta_int)^2
variance_all_cv <- apply(coeffs_all_cv, MARGIN = 2, var)
mse_all_cv <- bias2_all_cv + variance_all_cv

total_bias2_cv <- sum(bias2_all_cv)
total_variance_cv <- sum(variance_all_cv)
total_mse_cv <- total_bias2_cv + total_variance_cv

prob_true_model_cv <- mean(check_true_model_cv_all)
mean_nonzero_coeffs_cv <- mean(n_nonzero_coeffs_cv_all)

mean_lambda_min_cv <- mean(lambda_min_cv_all)

# Statistics for Lasso-OLS-CV
beta_int <- c(0, beta)
bias2_all_ols <- (colMeans(coeffs_all_ols) - beta_int)^2
variance_all_ols <- apply(coeffs_all_ols, MARGIN = 2, var)
mse_all_ols <- bias2_all_ols + variance_all_ols

total_bias2_ols <- sum(bias2_all_ols)
total_variance_ols <- sum(variance_all_ols)
total_mse_ols <- total_bias2_ols + total_variance_ols

prob_true_model_ols <- mean(check_true_model_ols_all)
mean_nonzero_coeffs_ols <- mean(n_nonzero_coeffs_ols_all)

mean_lambda_min_ols <- mean(lambda_min_ols_all)

# Statistics for BIC
beta_int <- c(0, beta)
bias2_all_bic <- (colMeans(coeffs_all_bic) - beta_int)^2
variance_all_bic <- apply(coeffs_all_bic, MARGIN = 2, var)
mse_all_bic <- bias2_all_bic + variance_all_bic

total_bias2_bic <- sum(bias2_all_bic)
total_variance_bic <- sum(variance_all_bic)
total_mse_bic <- total_bias2_bic + total_variance_bic

prob_true_model_bic <- mean(check_true_model_bic_all)
mean_nonzero_coeffs_bic <- mean(n_nonzero_coeffs_bic_all)

mean_lambda_min_bic <- mean(lambda_min_bic_all)
# (Bias, variance, MSE calculations, plotting, etc.)

# get MSE for OLS on only relevant variables
beta_int <- c(0, beta)
bias2_all_ols_true <- (colMeans(coeffs_ols_true) - beta_int)^2
variance_all_ols_true <- apply(coeffs_ols_true, MARGIN = 2, var)
mse_all_ols_true <- bias2_all_ols_true + variance_all_ols_true

total_bias2_ols_true <- sum(bias2_all_ols_true)
total_variance_ols_true <- sum(variance_all_ols_true)
total_mse_ols_true <- total_bias2_ols_true + total_variance_ols_true



# Summary of results
results_summary <- data.frame(MSE = c(total_mse_cv, total_mse_ols, total_mse_bic),
                              Prob_True_Model = c(prob_true_model_cv, prob_true_model_ols, prob_true_model_bic),
                              Mean_Nonzero_Coeffs = c(mean_nonzero_coeffs_cv, mean_nonzero_coeffs_ols, mean_nonzero_coeffs_bic),
                              Mean_Lambda_Min = c(mean_lambda_min_cv, mean_lambda_min_ols, mean_lambda_min_bic)
)  %>% 
  rbind(total_mse_bench)  %>% 
    rbind(total_mse_ols_true)
row.names(results_summary) <- c("Lasso-CV", "Lasso-OLS-CV", "Lasso BIC", "OLS-20", "OLS-5")
results_summary[4,2:4] <- c(0, k, NA)
results_summary[5,2:4] <- c(0, k_rel, NA)
print(results_summary)
#%%

getwd()
write.csv(results_summary, file = "Selection/data/selection.csv")
#write.csv(results_summary, file = file.path(getwd(), "comparison_cv_likelihood.csv"))
#write.csv(results_summary, file = file.path(getwd(), "comparison_cv_likelihood.csv"))

stargazer(results_summary, summary = FALSE, rownames = TRUE, title = "Comparison of Parameter Selection Procedures", out = "Selection/tables/selection.tex")

#%%
