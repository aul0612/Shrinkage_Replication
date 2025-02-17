rm(list = ls())

library(ggplot2)

# Define the Laplace prior (Lasso prior)
laplace_prior <- function(beta, lambda) {
  (lambda / 2) * exp(-lambda * abs(beta))
}

# Define the Gaussian prior (Ridge prior) for comparison
gaussian_prior <- function(beta, sigma) {
  (1 / sqrt(2 * pi * sigma^2)) * exp(-beta^2 / (2 * sigma^2))
}

beta_values <- seq(-4, 4, length.out = 500)

lambda <- 2  # Laplace parameter
sigma <- 1 / sqrt(2 * lambda)  
laplace_density <- laplace_prior(beta_values, lambda)
gaussian_density <- gaussian_prior(beta_values, sigma)

data <- data.frame(
  beta = rep(beta_values, 2),
  density = c(laplace_density, gaussian_density),
  prior = rep(c("Laplace (Lasso)", "Gaussian (Ridge)"), each = length(beta_values))
)

# Only data entries for Lasso
data_lasso <- subset(data,prior == "Laplace (Lasso)")

p <- ggplot(data_lasso, aes(x = beta, y = density)) +
  geom_line(size = 1) +
  labs(
    #title = "Bayesian Viewpoint on the Lasso Estimate",
    #subtitle = "Laplace Prior (Lasso)",
    x = expression(beta),
    y = "Density",
    color = "Prior"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),
    panel.grid.major.x = element_line(color = "gray", size = 0.5),
    panel.grid.major.y = element_line(color = "gray", size = 0.5), 
    panel.grid.minor.x = element_line(color = "lightgray", size = 0.25), 
    panel.grid.minor.y = element_line(color = "lightgray", size = 0.25)
  )

ggsave("Lasso/Output/bayesian_view_on_lasso.pdf",
       width = 16, height = 9, units = "cm",
       plot = p,
       device = cairo_pdf)