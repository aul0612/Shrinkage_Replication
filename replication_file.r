# Purpose: Replicate the simulations for the paper
#%%
### Check and install the required packages
packages <- c("MASS", "parallel", "Matrix", "progress", "caret", "ggpubr",
             "ggplot2", "dplyr", "glmnet", "mosaic", "tidyr", "gridExtra",
             "cv", "stargazer", "doMC", "foreach", "here", "lme4", "doParallel")
suppressMessages(install.packages(packages, dependencies = TRUE))
lapply(packages, suppressPackageStartupMessages(require), character.only = TRUE)

setwd(here()) 
#%%
### Section 2
# Replicate simulations for ridge regression
path <- "Ridge/Code"
files <- list.files(path = path)
for (file in files){
    path <- "Ridge/Code"
    print(paste("Running", file))
    source(file.path(getwd(), path, file))
}
rm(list = ls())
#%%
### Section 4
# Replicate simulations for Comparison of Ridge vs. Lasso
path <- "Comparison/Code"
files <- list.files(path = path)
for (file in files){
    path <- "Comparison/Code"
    if (file != "comparison_plots.r"){
        print(paste("Running", file))
        source(file.path(getwd(), path, file))
    }
}
source("Comparison/Code/comparison_plots.r")
rm(list = ls())
#%%
### Section 5
# Replicate simulations for tables in section 5.4
os <- .Platform$OS.type
if (os == "unix"){
    print("Running selection_unix.r")
    source("Selection/Code/selection_unix.r")
} else {
    print("Running selection_windows.r")
    source("Selection/Code/selection_windows.r")
}
#%%