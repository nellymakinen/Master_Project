#----------------------------------------------
#----------------------------------------------
# Script as pre-processing step for DIABLO
#
# Author: Nelly Mäkinen
# Date: 30th May, 2024
#----------------------------------------------
#----------------------------------------------

# Load necessary packages
library(stringr)
library(mixOmics) 
library(dplyr)

# Read necessary files
#bolomics <- read_tsv("C:/Users/nelly/Desktop/Skola/Masterprojekt/Results/Autoscaled_Data/autoscaled_hilic_pos_var0.05.tsv")
genomics <- readRDS("C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/pathway_specific_clr_dataset.rds")
genomics <- as.data.frame(genomics)
extradata <- read.csv("C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/240219_metadata.csv")
extradata <- as.data.frame(extradata)

# Create a new columns for CeGaT.ID using extradata
genomics$CeGaT.ID <- rownames(genomics)

# Set row names as sequential numbers
rownames(genomics) <- 1:nrow(genomics)

genomics <- genomics[c("CeGaT.ID", setdiff(names(genomics), "CeGaT.ID"))]

# Create sample.id column in genomics
genomics$sample.id <- NA
# Create Label column in genomics
genomics$Label <- NA

# Save the last digits to simplify sample.id
genomics$sample.id <- sub(".*?(\\d+)$", "\\1", genomics$sample.id)

# Fix "ID" in extradata
extradata$ID <- str_trim(extradata$ID)
extradata$CeGaT.ID. <- str_trim(extradata$CeGaT.ID.)

# Fix subgroup..group.and.timepoint. in extradata
n <- 0
for(i in 1:nrow(extradata)) {
  if (extradata$CeGaT.ID.[i] %in% genomics$CeGaT.ID)
  {
    n <- n+1
    genomics$sample.id[i] <- extradata$ID[i]
    genomics$Label[i] <- extradata$subgroup..group.and.timepoint.[i]
  }
}

genomics <- genomics[c("Label", setdiff(names(genomics), "Label"))]
genomics <- genomics[c("sample.id", setdiff(names(genomics), "sample.id"))]

# Extract the ids in bolomics
#bolomics$sample.id <- sub(".*_(\\d{2}_\\d{3})_.*", "\\1", bolomics$sample.id)
#bolomics$sample.id <- sub(".*fm_(\\d{2}_\\d{2}).*", "\\1", bolomics$sample.id)
#bolomics$sample.id <- sub(".*PTE_(\\d{2}_\\d{2}).*", "\\1", bolomics$sample.id)
#bolomics$sample.id <- sub(".*TBI_(\\d{2}_\\d{2}).*", "\\1", bolomics$sample.id)
#bolomics$sample.id <- sub(".*Fm_(\\d{2}_\\d{2}).*", "\\1", bolomics$sample.id)
#bolomics$sample.id

# For-loop to change format of Label in genomics
# Iterate over each row in 'genomics'
for(i in 1:nrow(genomics)) {
  # Check if 'Label' column contains "X_x"
  
  if("CTR_early" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "10d_ctr"
  }
  if("CTR_late" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "6_ctr"
  }
  if("PTE_early" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "10d_pte"
  }
  if("PTE_late" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "6_PTE"
  }
  if("TBI_early" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "10d_TBI"
  }
  if("TBI_late" %in% genomics$Label[i]) {
    # Update the label in 'genomics'
    genomics$Label[i] <- "6_TBI"
  }
}

genomics <- genomics %>%
  select(sample.id, Label, everything())

# Save the pre-processed files
write_tsv(genomics, "C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/clr_data_fixed_for_DIABLO_specific_pathways.tsv")
#write_tsv(bolomics, "C:/Users/nelly/Desktop/Skola/Masterprojekt/MetabolomicsData/manipulated_datasets/ForDIABLO/var_threshold_0.05/hilic_pos_data_fixed_for_DIABLO.tsv")


