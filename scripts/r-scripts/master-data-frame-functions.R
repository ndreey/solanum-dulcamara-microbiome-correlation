library(tidyverse)

# Author: André Bourbonnais
# email: andbou95@gmail.com
# Last update: 2025-12-17
#
# This script holds functions for the 02-master-data-frame.R script.
# To make the exploration and mantel tests easier.
#

set.seed(1337)

#### FUNCTIONS #####
# calculates the mean bray-curtis distances 
get_bray_curtis <- function(abund_matrix, dataset_name, meta_data) {
  
  message(paste("Processing", dataset_name, "..."))
  
  # Transpose and calculate Bray-Curtis
  asv <- abund_matrix %>% t()
  bray_dist <- vegdist(asv, method = "bray")
  
  # Convert to long format with population info
  df_bray <- as.matrix(bray_dist) %>%
    as.data.frame() %>%
    rownames_to_column("sample1") %>%
    pivot_longer(-sample1, names_to = "sample2", values_to = "bray") %>%
    left_join(meta_data %>% select(sampleID, pop.id), 
              by = c("sample1" = "sampleID")) %>%
    rename(pop1.id = pop.id) %>%
    left_join(meta_data %>% select(sampleID, pop.id), 
              by = c("sample2" = "sampleID")) %>%
    rename(pop2.id = pop.id) %>%
    mutate(pair = paste0(pop1.id, "-", pop2.id))
  
  # Calculate mean/median per pair
  df_bray_stats <- df_bray %>%
    group_by(pair) %>%
    summarise(
      mean_bray = mean(bray),
      median_bray = median(bray),
      n = n(),
      .groups = "drop"
    )
}

# Function to create distance matrix and run Mantel test on any of the columns
# of "df". For example, we can compare mean_bray.AMF vs mean_bray.ITS
run_mantel <- function(df, measure_y, measure_x) {
  
  y_dist <- df %>% 
    select(pop1.id, pop2.id, all_of(measure_y)) %>% 
    pivot_wider(names_from = pop2.id, values_from = all_of(measure_y)) %>% 
    column_to_rownames("pop1.id") %>% 
    as.dist()
  
  
  x_dist <- df %>% 
    select(pop1.id, pop2.id, all_of(measure_x)) %>% 
    pivot_wider(names_from = pop2.id, values_from = all_of(measure_x)) %>% 
    column_to_rownames("pop1.id") %>% 
    as.dist()
  
  mantel_res <- mantel(y_dist, x_dist)
  
  list(
    mantel = mantel_res,
    r = round(mantel_res$statistic, 4),
    p = round(mantel_res$signif, 4)
  )
}