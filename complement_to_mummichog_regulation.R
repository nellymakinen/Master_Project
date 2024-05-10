# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# Read the Data
metabolite_data <- read_tsv("your_data.tsv")

# Define your groups

group1 <- 1:5  
group2 <- 6:10  

# Function to determine regulation status
determine_regulation <- function(t_score, p_value, alpha = 0.05) {
  if (p_value < alpha) {
    if (t_score > 0) {
      return("up")
    } else if (t_score < 0) {
      return("down")
    }
  }
  return("not_regulated")
}

# Analyze each metabolite
results <- metabolite_data %>%
  mutate(
    t_score = map2_dbl(
      .x = select(., group1),
      .y = select(., group2),
      .f = ~ t.test(.x, .y)$statistic
    ),
    p_value = map2_dbl(
      .x = select(., group1),
      .y = select(., group2),
      .f = ~ t.test(.x, .y)$p.value
    ),
    regulation = pmap_chr(
      list(t_score, p_value),
      function(t, p) determine_regulation(t, p)
    )
  )

# Select relevant columns to write to the output
output_data <- results %>%
  select(Metabolite = 1, t_score, p_value, regulation)

# Write the results to a new TSV file
write_tsv(output_data, "metabolite_regulation_results_c18pos.tsv")

# Print out a message
cat("Analysis complete, results written to 'metabolite_regulation_results.tsv'\n")
