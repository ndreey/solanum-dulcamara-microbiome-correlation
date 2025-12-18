library(tidyverse)
library(vegan)
library(writexl)


# Author: André Bourbonnais
# email: andbou95@gmail.com
# Last update: 2025-12-18
#
# This scripts runs mantel tests on the columns vs Fst and stores it as .
# and excel file (.xlsx).
# 
# Note, although seed is set, re-running mantel tests will give different 
# p-values due to the permutations of the algorithm.
#


set.seed(1337)
setwd("../solanum-dulcamara-microbiome-correlation/")

source("scripts/r-scripts/master-data-frame-functions.R")

# Load data
master_frames <- readRDS("data/master_frames.Rds")

# Define taxonomic levels
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Define the microbiome types and their corresponding column names
microbiome_types <- list(
  "16S" = "mean_bray.16S",
  "ITS" = "mean_bray.ITS",
  "AMF" = "mean_bray.AMF",
  "GEO" = "km.diff"
)

# Initialize results dataframe
results <- data.frame(taxa = taxa_levels)

# Loop through each microbiome community
for (mb_name in names(microbiome_types)) {
  mb_col <- microbiome_types[[mb_name]]
  
  r_values <- c()
  p_values <- c()
  
  # Loop through each taxonomic level
  for (taxa in taxa_levels) {
    df <- master_frames[[taxa]]
    
    # Check if the microbiome column exists in the data
    if (mb_col %in% colnames(df)) {
      # Run Mantel test
      mantel_result <- run_mantel(df, mb_col, "fst")
      r_values <- c(r_values, round(mantel_result$r, 3))
      p_values <- c(p_values, round(mantel_result$p, 3))
    } else {
      # If column doesn't exist, add NA
      r_values <- c(r_values, NA)
      p_values <- c(p_values, NA)
    }
  }
  
  # Add columns to results
  results[[paste0("mantel.r_", mb_name)]] <- r_values
  results[[paste0("mantel.p_", mb_name)]] <- p_values
}

# View results
print(results)

# Write to Excel
write_xlsx(results, "mantel_test_results.xlsx")










