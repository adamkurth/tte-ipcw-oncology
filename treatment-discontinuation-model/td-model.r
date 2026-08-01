# # Target Trial Emulation (TTE) using Inverse Probability Clone Censoring Weighting (IPCW) and Regularized Regression of Flatiron database for advanced NSCLC patients treated with immunotherapy

# treatment-discontinuation-model
# When is the optimal time to discontinue ICI therapy?
################################# 
rm(list=ls()); set.seed(20260717)

library(ggplot2); library(dplyr) # visual
library(glmnet); library(splines); library(survival) # modeling

dir.path <- file.path("treatment-discontinuation-model")
fig.path <- file.path("figures", "td-model")

data <- readRDS(file.path(dir.path, "gp-data.rds")) # use same data for now
X <- data$X; X.long <- data$X.long; X.lc <- data$X.lc








