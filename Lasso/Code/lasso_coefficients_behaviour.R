rm(list = ls())

library(glmnet)
library(tidyr)

# Generate data
set.seed(12)
n <- 500    
p <- 10       
X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
beta <- c(runif(7, -5, 5), rep(0, p - 7)) 
y <- X %*% beta + rnorm(n,0,2)
lambdas <- seq(0, 5.5, 0.25) 

fit <- glmnet(X, y, alpha = 1, lambda = lambdas)  # alpha=1 for Lasso regression

lambda <- fit$lambda                      
coefficients <- as.matrix(fit$beta)       
coef_df <- data.frame(lambda = lambda, t(coefficients))  # Data frame
coef_df_long <- pivot_longer(coef_df, -lambda, 
                             names_to = "Variable", 
                             values_to = "Coefficient")

true_beta_df <- data.frame(
  Variable = paste0("V", 1:p),
  TrueBeta = beta
)

# Plot using ggplot2
p <- ggplot(coef_df_long, aes(x = lambda, y = Coefficient, color = Variable)) +
  geom_line(size = 0.7) +
  labs(
       x = expression(lambda),
       y = "Coefficients") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_discrete(name = "Predictors") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_x_continuous(breaks = seq(0, 6, 0.5))

ggsave("Lasso/Output/coefficients_behaviour.pdf",
       width = 16, height = 9, units = "cm",
       plot = p,
       device = cairo_pdf)
