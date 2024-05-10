##--------------------------##
## Code to implement GBMs   ##
##--------------------------##
## I don't know what you    ##
## need or want but I       ##
## included parts I got     ##
## from AI, thought it      ##
## explains more :)         ##
## Hope this is relevant!   ##
##--------------------------##

install.packages(c("tidyverse", "vegan", "ggplot2", "biomformat", "phyloseq", "igraph", "ggraph"))
library(tidyverse)
library(vegan)
library(ggplot2)
library(biomformat)
library(phyloseq)
library(igraph)
library(ggraph)

# Download gbms package
# This didn't work for me so I downloaded manually from:
# https://raeslab.org/software/gbms.html
download.file("https://raeslab.org/software/gbms.html", 
              destfile = "GBMs.zip")

# Unzip the GBM folder manually and place it outside the zipped file

# Check what the folder includes
list.files("C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/GBMs")

# Read GBMs data 
gbms_data <- read.csv("C:/Users/nelly/Desktop/Skola/Masterprojekt/metagenome/GBMs/GBM.inputfile.txt")

# Quick look at the data structure
head(gbms_data)

# Summary of the data
summary(gbms_data)

##-------------------------------------------------
# I don't know if you need these parts but include 
# them anyway hihi, just in case
# It's a bit of pre-processing, which I guess
# you don't need but I let it be left in the script
# so I don't accidentally remove an important step!
##-------------------------------------------------

# Load BIOM-format metagenomics data
biom_data <- read_biom("path/to/your/otu_table.biom")

# Convert to a phyloseq object, if you have taxonomy and sample data
ps <- phyloseq(otu_table(biom_data, taxa_are_rows = FALSE))

# Example of transforming counts to relative abundance
ps_relab <- transform_sample_counts(ps, function(x) x / sum(x))

# Optionally, filter rare taxa to reduce noise
ps_relab <- prune_species(speciesSums(ps_relab) > 0.001, ps_relab)

##-------------------------------------------
# Correlation analysis
##-------------------------------------------

# Example: Match GBM metabolites with metagenomic features
matched_data <- gbms_data %>%
  filter(Metabolite %in% taxa_names(ps_relab))

# Compute correlations
cor_results <- cor(as.matrix(otu_table(ps_relab)), as.matrix(matched_data[ , 3:ncol(matched_data)]))

# Convert to a data frame for easier manipulation and visualization
cor_df <- as.data.frame(as.table(cor_results))

# Filter to show strong correlations only
significant_correlations <- cor_df %>% 
  filter(abs(Freq) > 0.6)

##------------------------------------------
## Multivariate Analysis
##------------------------------------------

# Prepare a data frame of GBMs for RDA or CCA
gbms_env <- as.data.frame(matched_data[ , 3:ncol(matched_data)])

# Run RDA
rda_result <- rda(otu_table(ps_relab) ~ ., gbms_env)

# Plot the results
plot(rda_result)

##----------------------------------------
## Differential Abundance Analysis
##----------------------------------------
# Example using DESeq2 for differential abundance in pathways
library(DESeq2)

# Convert phyloseq to DESeq2
dds <- phyloseq_to_deseq2(ps_relab, ~ 1)
dds <- DESeq(dds)

# Results for pathways
res <- results(dds)
resLFC <- lfcShrink(dds, coef="condition_trt_vs_ctrl", type="normal")

# Identify significant pathways
sig_res <- resLFC %>% 
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  filter(padj < 0.05)

##---------------------------------------
## Visualization
##---------------------------------------

# Plot Significant Correlations

ggplot(significant_correlations, aes(x=Var1, y=Var2, fill=Freq)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title="Correlations between Microbes and GBMs", x="Microbial Taxa", y="GBM Metabolites")


# Network visualization

# Create a graph from correlations
graph_data <- significant_correlations %>% 
  filter(abs(Freq) > 0.6) %>%
  graph_from_data_frame()

# Plot the network
ggraph(graph_data, layout = 'fr') +
  geom_edge_link(aes(edge_alpha = Freq), show.legend = FALSE) +
  geom_node_point(color = 'lightblue', size = 5) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 5) +
  theme_void()




