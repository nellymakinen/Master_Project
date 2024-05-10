#----------------------------------------------------
#----------------------------------------------------
# R code for autoscaling each data set
#
# Author: Nelly Mäkinen
# Date: 30th May, 2024
#----------------------------------------------------
#----------------------------------------------------

# Install and load necessary packages 
#install.packages("readr")
library(readr)
#install.packages("mixOmics")
library(mixOmics)
#install.packages("pls")
library(pls)
# For reproducibility
set.seed(123)
#---------------------------------------------------

# Read in the metabolomics data from a TSV file
data <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/MetabolomicsData/normalized_data/pte_metab_hilic_neg_data_normalized.tsv")
#---------------------------------------------------
#
# PRE-PROCESSING (Autoscaling)
#
#---------------------------------------------------
# Exclude the first two non-numeric columns
numeric_data <- data[, -c(1, 2)]  

# Calculate variances and directly use them for filtering
variances <- apply(numeric_data, 2, var)

# Set variance threshold
var_threshold <- 0.025

# Filter based on the variance threshold
keep <- variances > var_threshold
numeric_data_clean <- numeric_data[, keep]

# Autoscale the numeric data
centered_data <- scale(numeric_data_clean)  # scale() function centers and scales

# Combine with sample.id and Label
final_data <- data.frame(sample.id = data$sample.id, Label = data$Label, centered_data)

# Write the autoscaled dataset to a TSV file
write_tsv(final_data, "C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/Autoscaled_Data_OG/autoscaled_hilic_neg_var0.025.tsv")