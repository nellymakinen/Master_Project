#----------------------------------------------------
#----------------------------------------------------
# R code for visualizing PCA for each data set
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
first_data <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/Autoscaled_Data_OG/autoscaled_hilic_neg_var0.025.tsv")

#first_data <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/MetabolomicsData/manipulated_datasets/ForDIABLO/var_threshold_0.05/c18_neg_data_fixed_for_DIABLO.tsv")
##----------------------------------------------------

n_samples <- nrow(first_data)
n_metabolites <- ncol(first_data)

# ----------------------------------------------------
#
# PRE-PROCESSING (PCA)
#
# ----------------------------------------------------
# Extract the id part for the C18_neg data set
extracted_ids <- sub(".*_(\\d{2}_\\d{3})_.*\\.raw", "\\1", first_data$sample.id)
extracted_ids <- sub(".*fm_(\\d{2}_\\d{2}).*\\.raw", "\\1", extracted_ids)
extracted_ids <- sub(".*PTE_(\\d{2}_\\d{2}).*\\.raw", "\\1", extracted_ids)
extracted_ids <- sub(".*TBI_(\\d{2}_\\d{2}).*\\.raw", "\\1", extracted_ids)
extracted_ids <- sub(".*Fm_(\\d{2}_\\d{2}).*\\.raw", "\\1", extracted_ids)
extracted_ids
# Extract the id part for c18_pos, hilic_pos, hilic neg
extracted_ids <- sub(".*_(\\d{2}_\\d{3})_.*", "\\1", first_data$sample.id)
extracted_ids <- sub(".*fm_(\\d{2}_\\d{2}).*", "\\1", extracted_ids)
extracted_ids <- sub(".*PTE_(\\d{2}_\\d{2}).*", "\\1", extracted_ids)
extracted_ids <- sub(".*TBI_(\\d{2}_\\d{2}).*", "\\1", extracted_ids)
extracted_ids <- sub(".*Fm_(\\d{2}_\\d{2}).*", "\\1", extracted_ids)
extracted_ids

# Create a data frame
data <- data.frame(sample.id = extracted_ids, label = rep(NA, length(extracted_ids)))
data[1:50,2] <- "early"
data[51:95,2] <- "late"
label <- as.factor(data[,2])

# Perform PCA coloring after time-point
pca.result_full <- pca(first_data[,-c(1,2)], ncomp = 2)  # ncomp is the number of principal components

# Direct mapping 
pca.result_full[["names"]][["sample"]] <- extracted_ids[as.numeric(pca.result_full[["names"]][["sample"]])]

plotIndiv(pca.result_full, 
          legend = TRUE,
          group = data$label,
          cex = 3,
          pch = 17,
          col = c("red4", "orange"),
          title = "Early and Late Labels - C18 Pos")
# -----------------------------------------------------
# PCA plots specifically for the early time point
# -----------------------------------------------------
subdata_low <- first_data[1:50,]
labels <- as.factor(subdata_low$Label)

pca.result_early <- pca(subdata_low[,-c(1,2)], ncomp = 2)

early_extracted_ids <- extracted_ids[1:50]
pca.result_early[["names"]][["sample"]] <- early_extracted_ids[as.numeric(pca.result_early[["names"]][["sample"]])]

plotIndiv(pca.result_early, 
          group = labels,
          cex = 3,
          pch = 17,
          legend = TRUE,
          col = c("red4", "orange", "purple"),
          title = "PCA Plot early")
# -----------------------------------------------------

# Extract the ids part c18_pos, hilic_pos
data$sample.id <- sub(".*_(\\d{2}_\\d{3})_.*", "\\1", data$sample.id)
data$sample.id <- sub(".*fm_(\\d{2}_\\d{2}).*", "\\1", data$sample.id)
data$sample.id <- sub(".*PTE_(\\d{2}_\\d{2}).*", "\\1", data$sample.id)
data$sample.id <- sub(".*TBI_(\\d{2}_\\d{2}).*", "\\1", data$sample.id)
data$sample.id <- sub(".*Fm_(\\d{2}_\\d{2}).*", "\\1", data$sample.id)

# -----------------------------------------------------
# PCA specifically for late time point
#------------------------------------------------------
subdata_high <- first_data[48:89,]
labels <- as.factor(subdata_high$Label)

pca.result_late <- pca(subdata_high[,-c(1,2)], ncomp = 2)

late_extracted_ids <- extracted_ids[48:89]
pca.result_late[["names"]][["sample"]] <- late_extracted_ids

plotIndiv(pca.result_late, 
          group = labels, 
          legend = TRUE,
          title = "PCA Plot with CTR, PTE & TBI")
# -----------------------------------------------------

# First, create the PCA plot without individual names
plotIndiv(pca.result_full, 
          group = labels_factor_full_y, 
          ind.names = FALSE,
          legend = TRUE,     
          title = 'PCA Plot all samples Hilic Pos')
          #col = c("green", "purple", "orange"))

##----------------------------------------------------
##----------------------------------------------------

# Perform PCA coloring after time-point
pca.result_full <- pca(first_data, ncomp = 4)  # ncomp is the number of principal components

# Extract the variance explained by each principal component
var_explained <- pca.result_full$sdev^2
var_explained_percent <- var_explained / sum(var_explained) * 100

# Create a scree plot
plot(var_explained_percent, xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     type = "b", pch = 19, main = "Scree Plot")


# Make sure 'labels' is a factor
labels_factor_full <- as.factor(data$Label)

# Plot Regular PCA

plotIndiv(pca.result_full, 
          group = labels_factor_full,
          ind.names = TRUE,  # Show sample numbers
          legend = TRUE,     # Add a legend if there are groups
          title = 'PCA Plot Hilic Pos',
          col = c("red", "blue", "green", "yellow", "purple", "orange"))

# Plot PCA results with colored groups for each label
plotIndiv(pca.result_full, 
          group = labels_factor_full, 
          ind.names = TRUE,  
          legend = TRUE,     
          title = 'PCA Plot all samples Hilic Pos',
          col = c("red", "blue", "green", "yellow", "purple", "orange"),
          pch = 20) 

pca.result_early <- pca(first_data[1:52,], ncomp = 5)
labels_factor_early <- as.factor(data$Label[1:52])

# Plot regular PCA early
plotIndiv(pca.result_early, 
          group = labels_factor_early,
          ind.names = TRUE,  # Show sample numbers
          legend = TRUE,     # Add a legend if there are groups
          title = 'PCA Plot Hilic Pos',
          col = c("yellow", "purple", "orange"))

# Plot PCA results with colored groups for each label
plotIndiv(pca.result_early, 
          group = labels_factor_early, 
          ind.names = TRUE,  
          legend = TRUE,     
          title = 'PCA Plot early samples Hilic Pos',
          col = c("yellow", "purple", "orange"),
          pch = 20) 

pca.result_late <- pca(first_data[53:97,], ncomp = 5)
labels_factor_late <- as.factor(data$Label[53:97])

# Plot regular PCA late
plotIndiv(pca.result_late, 
          group = labels_factor_late,
          ind.names = TRUE,  # Show sample numbers
          legend = TRUE,     # Add a legend if there are groups
          title = 'PCA Plot Hilic Pos',
          col = c("yellow", "purple", "orange"))

# Plot PCA results with 3 late labels

plotIndiv(pca.result_late, 
          group = labels_factor_late, 
          ind.names = TRUE,  
          legend = TRUE,     
          title = 'PCA Plot late samples Hilic Pos',
          col = c("yellow", "purple", "orange"),
          pch = 20) 

#---------------------------------------
#---------------------------------------
##--------------------------------------
##--------------------------------------

# Perform PCA coloring after year
# Extract the year part
extracted_years <- sub(".*PTE_(\\d{2})_.*", "\\1", data$sample.id)
extracted_years <- sub(".*fm_(\\d{2})_.*", "\\1", extracted_years)
extracted_years <- sub(".*TBI_(\\d{2})_.*", "\\1", extracted_years)
extracted_years <- sub(".*Fm_(\\d{2})_.*", "\\1", extracted_years)

# Print the extracted years
print(extracted_years)

# Make sure 'extracted_years' is a factor
labels_factor_full_y <- as.factor(extracted_years)

# First, create the PCA plot without individual names
plotIndiv(pca.result_full, 
          group = labels_factor_full_y, 
          ind.names = TRUE,  # Set to FALSE as we will add custom text
          legend = TRUE,     
          title = 'PCA Plot all samples Hilic Pos',
          col = c("green", "purple", "orange"))

# First, create the PCA plot without individual names
plotIndiv(pca.result_full, 
          group = labels_factor_full_y, 
          ind.names = FALSE,  # Set to FALSE as we will add custom text
          legend = TRUE,     
          title = 'PCA Plot all samples Hilic Pos',
          col = c("green", "purple", "orange"),
          pch = 20)

labels_factor_early_y <- as.factor(extracted_years[1:52])
# Plot PCA results with colored groups for each label
plotIndiv(pca.result_early, 
          group = labels_factor_early_y, 
          ind.names = TRUE,  
          legend = TRUE,     
          title = 'PCA Plot early samples Hilic Pos',
          col = c("green", "purple"),
          pch = 20) 

# Plot PCA results with 3 late labels

labels_factor_late_y <- as.factor(extracted_years[53:97])

plotIndiv(pca.result_late, 
          group = labels_factor_late_y, 
          ind.names = FALSE,  
          legend = TRUE,     
          title = 'PCA Plot late samples Hilic Pos',
          col = c("yellow", "purple", "orange"),
          pch = 20) 

#---------------------------------------

