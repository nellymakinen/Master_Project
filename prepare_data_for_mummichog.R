#------------------------------------------------
#------------------------------------------------
# Script which prepares a .txt file for mummichog
#
# Author: Nelly Mäkinen
# Date: 30th May, 2024
#------------------------------------------------
#------------------------------------------------

# Load necessary packages
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(mixOmics) 
library(broom)

# Read the TSV file
metabolites <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/Autoscaled_Data_OG/autoscaled_c18_neg_var0.025.tsv")
mz_rt_data <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/MetabolomicsData/annotations/c18_neg_metabo_annotations.tsv")

# Choose subset here 
# TBI vs PTE-----------
# EARLY --
# C18pos: 15:39
# C18neg: 15:52
# Hilicpos: 15:52
# Hilicneg: 14:50
# LATE --
# C18neg: 66:97
# C18pos: 52:82
# Hilicpos: 67:97
# Hilicneg: 65:95
# CTRL vs PTE-----------
# EARLY --
# C18pos: 1:26
# LATE --
# C18neg: 1:34
# CTRL vs TBI-----------
# EARLY --
# C18pos: c(1:14, 27:39)
# LATE --
# C18neg: c(1:14, 35:52)

metabolites <- metabolites[c(1:14, 35:52),]
metabolites <- t(metabolites)

# Select the relevant columns for mummichog from annotation files
# c(1,4,5) för C18pos, C18neg
# c(1,2,3) för hilicneg
# c(1,3,4) för hilicpos
mz_rt_data <- mz_rt_data[, c(1,4,5)]

# Convert to a data frame
your_data_df <- as.data.frame(metabolites)

# Convert row names to a new column
your_data_df$metabolite <- rownames(your_data_df)

# Reorder to make 'metabolite' the first column
your_data_df <- your_data_df[, c("metabolite", setdiff(names(your_data_df), "metabolite"))]

metabolites <- your_data_df[-c(1),]

filtered_mzrt <- mz_rt_data %>%
  filter(metabolite.id %in% metabolites$metabolite)

metabolites <- metabolites[-c(1), -c(1)]

# Create columns for t-score and p-value
filtered_mzrt$t_score <- NA_real_
filtered_mzrt$p_value <- NA_real_

#------------TBI vs PTE----------------
# EARLY 
# C18pos: 1:12 vs 13:25
# C18neg: 1:20 vs 21:38
# Hilicpos: 1:20 vs 21:38
# Hilicneg: 1:19 vs 20:37
# LATE 
# C18pos: 1:16 vs 17:32
# C18neg: 1:15 vs 16:31
# Hilicpos: 1:15 vs 16:31
# Hilicneg: 1:15 vs 16:31
#------------CTRL vs PTE---------------
# EARLY
# C18pos: 1:14 vs 15:26
# LATE
# C18neg: 1:14 vs 15:34
#------------CTRL vs TBI---------------
# EARLY
# C18pos: 1:14 vs 15:27
# LATE
# C18neg: 1:14 vs 15:32

# Loop through each row of metabolites to calculate t_scores and p_values
for (i in 1:nrow(metabolites)) {
  # Extract data for current row for both groups
  group1_data <- as.numeric(metabolites[i, 1:14]) # Adjust indices as necessary
  group2_data <- as.numeric(metabolites[i, 15:32]) # Adjust indices as necessary
  
  # Perform t-test
  t_test_result <- t.test(group1_data, group2_data, na.action = na.exclude)
  
  # Store the results in mz_rt_data
  filtered_mzrt$t_score[i] <- t_test_result$statistic
  filtered_mzrt$p_value[i] <- t_test_result$p.value
}

# Check the first few rows to verify
head(filtered_mzrt)

reordered_data <- filtered_mzrt %>%
  dplyr::select(Average.Mz, `Average.Rt(min)`, p_value, t_score, metabolite.id)

# Rename columns to match the exact names specified 
reordered_data <- rename(reordered_data,
                         mz = Average.Mz,
                         rt = `Average.Rt(min)`,
                         `p-value` = p_value,
                         `t-score` = t_score,
                         `metabolite.id` = metabolite.id)

# View the reordered dataframe
print(reordered_data)

# Convert from min to sec if necessary
reordered_data$rt <- 60*(reordered_data$rt)

# Write the results to a new file
write.table(reordered_data,"C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/mummichog/version2/late/c18neg/c18neg_late_w_ctrl_tbi.txt" , row.names = FALSE, quote = FALSE)

# Set the path to your original file and the output file
input_file_path <- "C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/mummichog/version2/late/c18neg/c18neg_late_w_ctrl_tbi.txt"
output_file_path <- "C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/mummichog/version2/late/c18neg/c18neg_late_w_ctrl_tbi_fixed.txt"

# Read the original file
lines <- readLines(input_file_path)

# Replace spaces with tabs
formatted_lines <- sapply(lines, function(x) gsub(" ", "\t", x))

# Write the formatted lines to a new file
writeLines(formatted_lines, output_file_path)

# -------------------------------------
# Adjust p-value
# -------------------------------------

# Adjust p-values for multiple testing
reordered_data$adjusted_p_value <- p.adjust(reordered_data$`p-value`, method = "BH", n = length(reordered_data$`p-value`))

# Plot histogram of adjusted p-values
hist(reordered_data$adjusted_p_value, breaks = 100, main = "Adjusted P-value Distribution", col= c("red4"))

#----------------------------------------------
#----------------------------------------------


