# Purpose: Replicate the simulations for the paper
# Check and install the required packages
packages <- c("MASS", "parallel", "Matrix", "progress", "caret",
             "ggplot2", "dplyr", "glmnet", "mosaic", "tidyr", "gridExtra",
             "cv", "stargazer", "doMC", "foreach", "here", "lme4", "doParallel")
for (package in packages){
    if (!require(package)){
        install.packages(package)
        library(package)
    }
}

setwd(file.path(here(), "Replication"))  # Works outside RStudio
# Replicate simulations for ridge regression
path <- "Ridge/Code"
files <- list.files(path = path)
for (file in files){
    print(paste("Running", file))
    source(file.path(getwd(), path, file))
}

# Replicate simulations for Comparison of Ridge vs. Lasso
files <- list.files(path = "Comparison/Code")
for (file in files){
    if file != "comparison_plots.r"{
        source(file)
    }
}
source("Comparison/Code/comparison_plots.r")

# Replicate simulations for parameter selection
os <- .Platform$OS.type
os
if (os == "unix"){
    source("Selection/Code/selection_unix.r")
} else {
    source("Selection/Code/selection_windows.r")
}