#------------------------------------
#------------------------------------
# DIABLO pipeline
#
# Author: Nelly Mäkinen
# Date: 30th May, 2024
#------------------------------------
#------------------------------------

# Import necessary packages
library(readr)
library(mixOmics) 

# Read your files for the two-omics analysis
bolomics <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/Autoscaled_Data_OG/autoscaled_c18_pos_var0.025.tsv")
genomics <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/clr_data_fixed_for_DIABLO_CORRECT.tsv")

# For reproducibility
set.seed(123) 

# ------------------------------------
# PRE-PROCESSING
# ------------------------------------

# Reorder genomics data after "Label"
genomics <- genomics[order(genomics$Label), ]
# Remove the unnecessary column (CeGat.ID)
genomics <- genomics[, -c(3)]

# ---------------------------------------------------
# Choose subdata here!
# ---------------------------------------------------

# BOLOMICS
# C18 Neg
# early (1:52), late (53:97), early_w/o_ctrl (15:52), late_w/o (66:97)

# C18 Pos
# early (1:39), late (40:82), early_w/o_ctrl (15:39), late_w/o (52:82)

# Hilic Pos
# early (1:52), late (53:97), early_w/o_ctrl (15:52), late_w/o (67:97)

# Hilic Neg
# early (1:50), late (51:95), early_w/o_ctrl (14:50), late_w/o (65:95)
# ---
# GENOMICS
# early (1:45), late (46:91), early_w/o (14:45), late_w/o (60:91)

# ---------------------------------------------------
X2 <- bolomics[15:39,]
X1 <- genomics[14:45,]  

X <- list(metagenomics = X1, metabolomics = X2)

# Check dimensions of the data
lapply(X, dim) 

#----------------------------------------------------

# Run PCA for Metabolomics
pca.metabolomics <- pca(X2[,-c(1,2)], ncomp = 2)
# Run PCA for Metagenomics
pca.metagenomics <- pca(X1[,-c(1,2)], ncomp = 2)

labelb <- as.factor(X2$Label) 
labelg <- as.factor(X1$Label)

# Extract sample.id from metabolomics to match metagenomics sample.id
extracted_id <- sub(".*_(\\d{2}_\\d{3})_.*", "\\1", X2$sample.id)
extracted_id <- sub(".*fm_(\\d{2}_\\d{2}).*", "\\1", extracted_id)
extracted_id <- sub(".*PTE_(\\d{2}_\\d{2}).*", "\\1", extracted_id)
extracted_id <- sub(".*TBI_(\\d{2}_\\d{2}).*", "\\1", extracted_id)
extracted_id <- sub(".*Fm_(\\d{2}_\\d{2}).*", "\\1", extracted_id)

# Rename sample-names in pca results for metabolomics for easier interpretation of plot
pca.metabolomics[["names"]][["sample"]] <- extracted_id[as.numeric(pca.metabolomics[["names"]][["sample"]])]
# Rename sample-names for metabolomics
X$metabolomics$sample.id <- extracted_id

# Plot PCA for Metabolomics
plotIndiv(pca.metabolomics, 
          group = labelb, 
          title = 'PCA Metabolomics',
          legend = TRUE)

# Plot PCA for Metagenomics
plotIndiv(pca.metagenomics, 
          group = labelg, 
          title = 'PCA Metagenomics',
          legend = TRUE)

# -----------------------------------------------

# Check for different number of rows in the two data sets
sapply(X, nrow)

# Extract sample IDs from each block (X1, X2)
sample_ids_list <- lapply(X, function(df) df$sample.id)

# Find common sample IDs across both blocks
common_sample_ids <- Reduce(intersect, sample_ids_list)

# Subset each block to include only the common samples
X_aligned <- lapply(X, function(df) {
  df_filtered <- df[df$sample.id %in% common_sample_ids, ]
  df_ordered <- df_filtered[match(common_sample_ids, df_filtered$sample.id), ]
  return(df_ordered)
})

# Save the aligned sample id's from the two blocks
aligned_sample_ids <- X_aligned[[1]]$sample.id

# Create a response vector that matches 
Y_aligned <- as.factor(X_aligned$metabolomics$Label)
summary(Y_aligned)

# Remove the first two columns from each block (sample.id, Label)
X_aligned_numeric <- lapply(X_aligned, function(df) df[, -c(1, 2)])

colnames(X_aligned_numeric$metagenomics) <- sub(".*\\|(t__[^|]+)$", "\\1", colnames(X_aligned_numeric$metagenomics))

# Ensure that the number of rows in each block matches the length of Y_aligned
if(any(sapply(X_aligned_numeric, nrow) != length(Y_aligned))) {
  stop("Mismatch in number of samples between X blocks and Y.")
}

# -----------------------------------------------
#         INITIAL ANALYSIS
#      Pairwise PLS Comparison
# -----------------------------------------------

# Select arbitrary values of features to keep
list.keepX = c(10, 80) 
list.keepY = c(10, 80)

# Double-check for no zero-variance metabolites for X dataset
constant_vars_X <- apply(X_aligned_numeric[["metabolomics"]], 2, var) == 0
X_filtered <- X_aligned_numeric[["metabolomics"]][, !constant_vars_X]

# Double-check for no zero-variance metabolites for Y dataset
constant_vars_Y <- apply(X_aligned_numeric[["metagenomics"]], 2, var) == 0
Y_filtered <- X_aligned_numeric[["metagenomics"]][, !constant_vars_Y]

set_filtered <- list(metabolomics = X_filtered, metagenomics = Y_filtered)

# Run PLS
pls <- spls(X_filtered, Y_filtered, keepX = list.keepX, keepY = list.keepY)

# Plot features of PLS
plotVar(pls, cutoff = 0.5, title = "bolomics vs genomics", 
        legend = c("bolomics", "genomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# Calculate correlation of bolomics and genomics
cor <- cor(pls$variates$X, pls$variates$Y)

# Create a square matrix filled with the correlation value
design = matrix(cor[1,1], ncol = length(set_filtered), nrow = length(set_filtered), 
                dimnames = list(names(set_filtered), names(set_filtered)))

# Set diagonal to 0 to indicate no correlation within each block
diag(design) <- 0 
design

#------------------------------------------------------
#
#           DIABLO
#
#------------------------------------------------------

# For reproducibility
set.seed(123)
# Form the basic DIABLO model
basic.diablo.model = block.splsda(X = set_filtered, 
                                  Y = Y_aligned, ncomp = 6, 
                                  design = design)

# -----------------------------------------------------
# Scatterplot of explained variance of each component
# -----------------------------------------------------

# Extract the values
y_values <- basic.diablo.model[["prop_expl_var"]][["metabolomics"]]

# Create x-axis values
x_values <- seq_along(y_values)

# Create the scatterplot
plot(x_values, y_values,
     xlab = "Component",
     ylab = "Proportion of Explained Variance",
     main = "Explained Variance by components in Metabolomics",
     pch = 19, col = "blue")

# ----------------------------------------------------
# Tuning the number of components
# ----------------------------------------------------

set.seed(55)

# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'loo', nrepeat = 10) 

# plot output of tuning
plot(perf.diablo, col = c("red4", "hotpink", "purple")) 

# Show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 
perf.diablo$WeightedVote.error.rate

# --------------------------------------------------
# Tuning the number of features
# This section is manually rerun until both "ncomp" 
# and "test.keepX" are optimized (as low but as accurate 
# as possible)
# --------------------------------------------------

# Set grid of values for each component to test
test.keepX = list (metabolomics = seq(1, 150, by = 20), 
                   metagenomics = seq(1, 200, by = 20))

# For paralellization
BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores()-1)

# For reproducibility
set.seed(97)

# Run the feature selection tuning
tune.EPI = tune.block.splsda(X = set_filtered, Y = Y_aligned, ncomp = 4, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 3, nrepeat = 10,
                              dist = "max.dist", BPPARAM = BPPARAM)

# Analyse error-rates to find optimal number of features to keep
error = tune.EPI$error.rate

# Set the optimal values of features to keep
list.keepX = tune.EPI$choice.keepX 
list.keepX

# ------------------------------------------------
# FINAL MODEL
# ------------------------------------------------

# Form optimized DIABLO model
final.diablo.model = block.splsda(X = set_filtered, Y = Y_aligned, ncomp = 4, 
                                  keepX = list.keepX, design = design)

# Check design matrix for final model
final.diablo.model$design 

# ------------------------------------------------
# PERFORMANCE OF THE FINAL MODEL
# ------------------------------------------------

# For reproducibility
set.seed(65)

# run component number tuning with repeated CV
perf.diablo.tuned = perf(final.diablo.model, validation = 'loo', 
                         nrepeat = 10) 

# Plot output of tuning
plot(perf.diablo.tuned) 

# Show the optimal choice for ncomp for each dist metric
perf.diablo.tuned$choice.ncomp$WeightedVote

#perf.diablo.tuned$MajorityVote.error.rate
perf.diablo.tuned$WeightedVote.error.rate

# ----------------------------------------------
# GENERAL PLOTS 
# ----------------------------------------------

plotDiablo(final.diablo.model, ncomp = 2, col = c("orange", "red4"))

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')

plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'Arrow plot')

# ---------------------------------------------
# VARIABLE PLOTS
# ---------------------------------------------

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

circosPlot(final.diablo.model, comp = 1:2, cutoff = 0.8, line = TRUE, 
           color.blocks= c('yellow', 'red4'), 
           color.cor = c("hotpink", "purple3"), linkWidth = c(1,3),
           size.labels = 1.0, 
           size.variables = 0.7)

# -------------------------------------------
# KEY FEATURE PLOT
# -------------------------------------------

plotLoadings(final.diablo.model, comp = 1, contrib = 'max', 
             method = 'median', legend.color = c("hotpink", "purple3"))



