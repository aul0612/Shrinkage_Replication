## data generation and metrics
#library(parallel)

#files <- list("Comparison/control_k_large.r", "Comparison/true_k_large.r", "Comparison/lasso_better.r", "Comparison/ridge_better.r")
# run each file in parallel
#cl <- makeCluster(4)
#parLapply(cl, files, source)
# stop the cluster
#stopCluster(cl)
#print("done")

## 
library(tidyverse)
library(ggplot2)


# read data
files <- list.files(path = "Comparison/data", pattern = "csv", full.names = TRUE)
data <- lapply(files, read.csv)
names(data) <- gsub("Comparison/data/", "", gsub(".csv", "", files))
names(data)

# lambda range for each simulation
lambda_values <- seq(0, 10, length.out = 100)

# adding lambda values to each dataframe in the list
data <- map(data, ~mutate(.x, lambda = lambda_values))


# transform data to long format in one dataframe with column lambda and the file name
data_long <- map(data, ~gather(.x, key = metric, value = value, -lambda))  %>% 
    bind_rows(.id = "file") %>% 
    as_tibble()


# data_long <- map_dfr(data, ~gather(.x, key = metric, value = value, -lambda))

# drop columns containing the coefficients
data_long <- data_long %>% 
  filter(!grepl("X", metric))


data_long_k20 <- data_long %>% 
  filter(file %in% c("lasso_better", "ridge_better"))
data_long_k100 <- data_long %>% 
  filter(file %in% c("control_k_large", "true_k_large"))

# Create the plot with two y-axes
# all metrics with "lasso" in dashed lines, all with "ridge" in solid lines
# mse in black, bias in red, variance in blue, mean-nonzero-coefficients in gray
# one plot for each file_name

##### plot the data for k=20
# creating a labeller
my_labeller <- c("ridge_better" = "(b)", "lasso_better" = "(a)")

plots_k20 <- data_long_k20 %>%
  ggplot(aes(x = lambda, y = value)) +
  # Layer for general metrics (mse, bias, variance)
  geom_line(aes(
    color = case_when(
      grepl("mse", metric) ~ "MSE",
      grepl("bias", metric) ~ "sq. Bias",
      grepl("variance", metric) ~ "Var",
      TRUE ~ "other"
    ), 
    linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
  )) +
  # Layer for 'mean-nonzero-coefficients' (on the secondary axis)
  geom_line(data = data_long_k20 %>% filter(grepl("mean_nonzero_coeffs", metric)),
            aes(
              x = lambda,
              y = value,
              color = "Nonzero Coefficients",  # Explicitly set color for legend
              linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
            ),  # Customize size if needed
            show.legend = TRUE) +  # Ensure it is included in the legend
  scale_color_manual(
    values = c(
      "MSE" = "black",
      "sq. Bias" = "red",
      "Var" = "blue",
      "Nonzero Coefficients" = "gray"  # Set color for Nonzero Coefficients in the legend
    )
  ) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("Lasso", "Ridge"), name = "Model Type") +
  facet_wrap(~file, scales = "free", labeller = labeller(file = my_labeller)) +  # Allow free scales for each facet
  labs(
    title = "",
    x = expression(lambda),
    y = "Metric Value",
    linetype = "Model Type",
    color = "Metric"
  ) +
  scale_y_continuous(
    limits = c(0, 7.5),  # Set limits for clarity
    sec.axis = sec_axis(~ ., name = "Mean Nonzero Coefficients")  # Adjust if needed
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 4.4)) +  # Set limits for clarity
  theme(
    strip.text = element_text(size = 10),   # Adjust strip text size for clarity
    #strip.background = element_blank(),  # Remove strip background for clarity
    panel.spacing = unit(1, "lines"),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 10))) 

# Display the plot
print(plots_k20)

ggsave("Comparison/plots/plotsk20.pdf",width = 16, height = 8, units = "cm")


##### plot the data for k=100
# creating a labeller
my_labeller <- c("control_k_large" = "(a)", "true_k_large" = "(b)")

plots_k100 <- data_long_k100 %>%
    # First filter out the nonzero_coeffs_lasso metric
    filter(!grepl("nonzero_coeffs_lasso", metric)) %>%
    ggplot(aes(x = lambda, y = value)) +
    # Layer for general metrics (mse, bias, variance)
    geom_line(aes(
        color = case_when(
            grepl("mse", metric) ~ "MSE",
            grepl("bias", metric) ~ "sq. Bias",
            grepl("variance", metric) ~ "Var",
            TRUE ~ "other"
        ), 
        linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
    )) +
    # Layer for 'mean-nonzero-coefficients' (on the secondary axis)
    geom_line(data = data_long_k100 %>% filter(grepl("mean_nonzero_coeffs", metric)),
                        aes(
                            x = lambda,
                            y = value*4,
                            color = "Nonzero Coefficients",  # Explicitly set color for legend
                            linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
                        ),  # Customize size if needed
                        show.legend = TRUE) +  # Ensure it is included in the legend
    scale_color_manual(
        values = c(
            "MSE" = "black",
            "sq. Bias" = "red",
            "Var" = "blue",
            "Nonzero Coefficients" = "gray50"  # Set color for Nonzero Coefficients in the legend
        )
    ) +
    scale_linetype_manual(
        values = c("solid", "dashed"),
        labels = c("Lasso", "Ridge"),
        name = "Model Type"
    ) +
    facet_wrap(~file, scales = "free", labeller = labeller(file = my_labeller)) +  # Allow free scales for each facet
    labs(
        title = "",
        x = expression(lambda),
        y = "Metric Value",
        color = "Metric",
        sec.axis = sec_axis(~ ., name = "Mean Nonzero Coefficients")  # Adjust if needed
    ) +
    scale_y_continuous(
        limits = c(0, 400),  # Set limits for clarity
        sec.axis = sec_axis(~ ./4, name = "Mean Nonzero Coefficients")  # Adjust if needed
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 10)) +  # Set limits for clarity
    theme(
        strip.text = element_text(size = 10),   # Adjust strip text size for clarity
        #strip.background = element_blank(),  # Remove strip background for clarity
        panel.spacing = unit(1, "lines"),
        strip.placement = "outside",
        axis.title.x = element_text(margin = margin(t = 10)))

# Display the plot
print(plots_k100)

#ggsave("Comparison/plots/plotsk100.pdf",width = 16, height = 8, units = "cm")




###### nicer plot with two seperate plots and a shared legend


# Create ggplot object for the "control_k_large" file
plot_control_k_large <- data_long_k100 %>%
    filter(file == "control_k_large", !grepl("nonzero_coeffs_lasso", metric)) %>%
    ggplot(aes(x = lambda, y = value)) +
    geom_line(aes(
        color = case_when(
            grepl("mse", metric) ~ "MSE",
            grepl("bias", metric) ~ "sq. Bias",
            grepl("variance", metric) ~ "Var",
            TRUE ~ "other"
        ), 
        linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
    )) +
    geom_line(data = data_long_k100 %>% 
                  filter(file == "control_k_large", grepl("mean_nonzero_coeffs", metric)),
              aes(
                  x = lambda,
                  y = value ,
                  color = "Nonzero Coefficients",  
                  linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
              ), 
              show.legend = TRUE) +
    scale_color_manual(
        values = c(
            "MSE" = "black",
            "sq. Bias" = "red",
            "Var" = "blue",
            "Nonzero Coefficients" = "gray50"
        )
    ) +
    scale_linetype_manual(
        values = c("solid", "dashed"),
        labels = c("Lasso", "Ridge"),
        name = "Model Type"
    ) +
    labs(
        title = "",
        x = expression(lambda),
        y = "Metric Value",
        color = "Metric",
        sec.axis = sec_axis(~ ., name = "")
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        sec.axis = sec_axis(~ ., name = "")
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 5)) +
    theme(
        #legend.position = "none",
        strip.text = element_text(size = 10),
        panel.spacing = unit(1, "lines"),
        axis.title.x = element_text(margin = margin(t = 10))
    )

# Create ggplot object for the "true_k_large" file
plot_true_k_large <- data_long_k100 %>%
    filter(file == "true_k_large", !grepl("nonzero_coeffs_lasso", metric)) %>%
    ggplot(aes(x = lambda, y = value)) +
    geom_line(aes(
        color = case_when(
            grepl("mse", metric) ~ "MSE",
            grepl("bias", metric) ~ "sq. Bias",
            grepl("variance", metric) ~ "Var",
            TRUE ~ "other"
        ), 
        linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
    )) +
    geom_line(data = data_long_k100 %>% 
                  filter(file == "true_k_large", grepl("mean_nonzero_coeffs", metric)),
              aes(
                  x = lambda,
                  y = value * 4,
                  color = "Nonzero Coefficients",  
                  linetype = ifelse(grepl("lasso", metric), "dashed", "solid")
              ), 
              show.legend = TRUE) +
    scale_color_manual(
        values = c(
            "MSE" = "black",
            "sq. Bias" = "red",
            "Var" = "blue",
            "Nonzero Coefficients" = "gray50"
        )
    ) +
    scale_linetype_manual(
        values = c("solid", "dashed"),
        labels = c("Lasso", "Ridge"),
        name = "Model Type"
    ) +
    labs(
        title = "",
        x = expression(lambda),
        y = "",
        color = "Metric",
        sec.axis = sec_axis(~ ., name = "Mean Nonzero Coefficients")
    ) +
    scale_y_continuous(
        limits = c(200, 400),
        sec.axis = sec_axis(~ ./4, name = "Mean Nonzero Coefficients")
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 5)) +
    theme(
        strip.text = element_text(size = 10),
        panel.spacing = unit(1, "lines"),
        axis.title.x = element_text(margin = margin(t = 10))
    )

library(gridExtra)
library(grid)
library(ggpubr)  # For get_legend function

# Function to extract legend
get_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]  # Extract the first item
    legend
}

# Extract legend from the first plot
shared_legend <- get_legend(plot_control_k_large)

# Remove legends from both plots
plot_control_k_large_no_legend <- plot_control_k_large +
    theme(legend.position = "none", )

plot_true_k_large_no_legend <- plot_true_k_large +
    theme(legend.position = "none")
title_control_k_large <- textGrob("(a)", gp = gpar(fontsize = 12))
title_true_k_large <- textGrob("(b)", gp = gpar(fontsize = 12))    

# Combine plots with the shared legend
plots_k100 <- grid.arrange(
    arrangeGrob(
        title_control_k_large,
        title_true_k_large, 
        plot_control_k_large_no_legend, 
        plot_true_k_large_no_legend, 
        nrow = 2,
        heights = c(1, 10), 
        widths = c(1, 1)  # Ensure equal widths
    ),
    shared_legend,  # Add the single legend as a grob
    ncol = 2,
    widths = c(10, 3.2)  # Allocate space for legend
)

ggsave("Comparison/plots/plotsk100.pdf",plots_k100, width = 16, height = 8, units = "cm")

