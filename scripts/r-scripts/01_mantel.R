library(tidyverse)
library(vegan)
library(scales)
library(ggpubr)
library(ggridges)
library(patchwork)
library(writexl)


# Author: André Bourbonnais
# email: andbou95@gmail.com
# Last update: 2025-12-17
#
# This scripts plots and summarizes most of the data to make comparisons.
# Showing distributions of the Bray-Curtis measurments, and plotting 
# results together.
#


#### FUNCTIONS #####
# FUNCTION TO PROCESS EACH DATASET
process_bray_curtis <- function(abund_matrix, dataset_name, meta_data) {
  
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
  
  # Ridge plot
  p_ridge <- df_bray %>%
    ggplot(aes(x = bray, y = reorder(pair, bray, mean), fill = after_stat(x))) +
    geom_density_ridges_gradient(alpha = 0.6, scale = 2) +
    geom_point(data = df_bray_stats, aes(x = mean_bray, y = pair), 
               color = "red", size = 2, alpha = 0.7, inherit.aes = FALSE) +
    geom_point(data = df_bray_stats, aes(x = median_bray, y = pair), 
               color = "blue", size = 2, alpha = 0.7, inherit.aes = FALSE) +
    scale_fill_viridis_c(guide = "none", option = "C") +
    labs(
      y = "Population Pair", 
      x = "Mean Bray-Curtis Distance",
      subtitle = "Red = mean, Blue = median"
    ) +
    theme_pubr(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))
  
  list(stats = df_bray_stats, ridge = p_ridge)
}

# CREATE FST VS BRAY PLOTS
create_fst_plots <- function(df, bray_col, dataset_name, mantel_res) {
  
  # Excluding self-comparisons
  p_no_self <- df %>%
    filter(pop1.id < pop2.id) %>% 
    ggplot(aes(x = .data[[bray_col]], y = fst)) +
    geom_smooth(method = lm, se = FALSE, color = "red", 
                linewidth = 0.9, linetype = "dashed") +
    geom_point(size = 2, alpha = 0.6) +
    labs(
      x = "Mean Bray-Curtis Distance", 
      y = "FST",
      subtitle = paste0(dataset_name,"\nMantel r = ", mantel_res$r, ", p = ", mantel_res$p)
    ) +
    theme_pubr(base_size = 12)
  
  # Including self-comparisons
  p_with_self <- df %>%
    filter(pop1.id <= pop2.id) %>%
    ggplot(aes(x = .data[[bray_col]], y = fst)) +
    geom_smooth(method = lm, se = FALSE, color = "red", 
                linewidth = 0.9, linetype = "dashed") +
    geom_point(size = 2, alpha = 0.6) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    labs(
      x = "Mean Bray-Curtis Distance", 
      y = "FST"
      ) +
    theme_pubr(base_size = 12)
  
  list(no_self = p_no_self, with_self = p_with_self)
}

# Function to create distance matrix and run Mantel test
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

##### PARAMETERS #####
set.seed(1337)
OUTDIR_RIDGE_AMF <- "plots/mantel-overview/ridge-combo-AMF"
OUTDIR_RIDGE_ITS <- "plots/mantel-overview/ridge-combo-ITS"
OUTDIR_RIDGE_16S <- "plots/mantel-overview/ridge-combo-16S"
OUTDIR_MANTEL_TOGETHER <- "plots/mantel-overview/mantel-together"
OUTDIR_MANTEL_AMF <- "plots/mantel-overview/mantel-AMF"
OUTDIR_MANTEL_ITS <- "plots/mantel-overview/mantel-ITS"
OUTDIR_MANTEL_16S <- "plots/mantel-overview/mantel-16S"

# Create output directory once
dir.create(OUTDIR_RIDGE_AMF, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_RIDGE_ITS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_RIDGE_16S, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_MANTEL_TOGETHER, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_MANTEL_AMF, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_MANTEL_ITS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_MANTEL_16S, showWarnings = FALSE, recursive = TRUE)


#### Loading Data  #####
# Load meta data
meta_dulcamara <- read_delim("doc/meta-microeco.tsv", delim = "\t", 
                             col_names = T) %>% 
  filter(spp == "Solanum_dulcamara")

meta_AMF <- read_delim("doc/metadata-AMF.tsv", delim = "\t", 
                       col_names = T) %>% 
  filter(spp == "Solanum_dulcamara") %>% 
  rename(pop.id = population)

# Load Hudson Fst data
df_piawka <- read_delim("04-piawka/final.merged.tsv",
                        delim = "\t", col_names = T)

# Load meta mapping pop -> short ids, habitat, humidity
pop_map <- read_delim("doc/uniq_pop_meta.tsv", delim = "\t", col_names = T)


#### Prepare data #####
# FST data
df_fst <- df_piawka %>% 
  select(pop1, pop2, metric, value) %>% 
  filter(!metric %in% c("pi", "Dxy")) %>% 
  rename(fst = value) %>% 
  select(pop1, pop2, fst) %>% 
  left_join(pop_map %>% select(pop, pop.id), by = c("pop1" = "pop")) %>%
  rename(pop1.id = pop.id) %>%
  left_join(pop_map %>% select(pop, pop.id), by = c("pop2" = "pop")) %>%
  rename(pop2.id = pop.id) %>% 
  mutate(pair = paste0(pop1.id, "-", pop2.id)) %>% 
  select(pair, pop1.id, pop2.id, fst)

# Make symmetric matrix of fst
df_sym <- df_fst %>% 
  bind_rows(
    df_fst %>% rename(pop1.id = pop2.id, pop2.id = pop1.id)
  ) %>% 
  bind_rows(
    tibble(
      pop1.id = unique(c(df_fst$pop1.id, df_fst$pop2.id)),
      pop2.id = unique(c(df_fst$pop1.id, df_fst$pop2.id)),
      fst = 0
    )
  ) %>%
  mutate(pair = paste0(pop1.id, "-", pop2.id))


# Initialize empty list to store results
mantel_summary <- tibble()

# Counter for taxa iteration tracking
counter <- 0

######## LOOP THROUGH TAXONOMIC LEVELS ########
for (TAXA in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  
  
  message(paste("\n========== Processing", TAXA, "=========="))
  
  counter <- counter + 1
  # Load abundance data
  abund_16S <- read_delim(paste0("finData/16S_abundance/", TAXA, "_abund.csv"),
                          delim = ",", col_names = T) %>% 
    rename("taxon" = "...1") %>% 
    filter(taxon != "k__Archaea") %>% 
    # Filter to only retain S. dulcamara samples
    column_to_rownames("taxon") %>%
    select(all_of(meta_dulcamara$sampleID))
  
  abund_ITS <- read_delim(paste0("finData/ITS_abundance/", TAXA, "_abund.csv"),
                          delim = ",", col_names = T) %>% 
    rename("taxon" = "...1") %>% 
    # Filter to only retain S. dulcamara samples
    column_to_rownames("taxon") %>%
    select(all_of(meta_dulcamara$sampleID))
  
  # AMF data is in different format, wrangle it to match
  abund_AMF <- read_delim(paste0("finData/AMF_abundance/", TAXA, "_abund.tsv"), 
                          delim = "\t", col_names = T, skip = 1) %>% 
    rename("taxon.tmp" = "#OTU ID") %>%
    # Split taxonomy string by semicolon
    separate(taxon.tmp, 
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";", 
             fill = "right",  # Fill missing columns with NA
             remove = FALSE) %>%  # Keep original taxon column
    # Replace empty strings with NA
    mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                  ~ifelse(. == "", NA, .))) %>%
    # Replace NA with appropriate prefix
    mutate(
      Kingdom = ifelse(is.na(Kingdom), "k__", Kingdom),
      Phylum = ifelse(is.na(Phylum), "p__", Phylum),
      Class = ifelse(is.na(Class), "c__", Class),
      Order = ifelse(is.na(Order), "o__", Order),
      Family = ifelse(is.na(Family), "f__", Family),
      Genus = ifelse(is.na(Genus), "g__", Genus),
      Species = ifelse(is.na(Species), "s__", Species)
    ) %>%
    # Create new standardized taxonomy string
    mutate(taxon = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "|")) %>%
    # Drop the parsed columns, keep only taxonomy and sample data
    select(taxon, starts_with("D")) %>% 
    # Filter to only retain S. dulcamara samples
    column_to_rownames("taxon") %>%
    select(all_of(meta_AMF$sampleID))


  #### PROCESS ALL DATASETS ####
  results_16S <- process_bray_curtis(abund_16S, "16S", meta_dulcamara)
  results_ITS <- process_bray_curtis(abund_ITS, "ITS", meta_dulcamara)
  results_AMF <- process_bray_curtis(abund_AMF, "AMF", meta_AMF)
  
  #### MERGE DATA ####
  df_mantel <- df_sym %>% 
    inner_join(results_AMF$stats, by = "pair") %>% 
    rename(mean_bray.AMF = mean_bray) %>% 
    inner_join(results_ITS$stats, by = "pair") %>% 
    rename(mean_bray.ITS = mean_bray) %>% 
    inner_join(results_16S$stats, by = "pair") %>% 
    rename(mean_bray.16S = mean_bray) %>% 
    select(pop1.id, pop2.id, fst, starts_with("mean_bray")) %>% 
    arrange(pop1.id, pop2.id)
  
  #### CREATE DISTANCE MATRICES AND RUN MANTEL TESTS ####
  # Run Mantel tests
  mantel_AMF <- run_mantel(df_mantel, "mean_bray.AMF", "fst")
  mantel_ITS <- run_mantel(df_mantel, "mean_bray.ITS", "fst")
  mantel_16S <- run_mantel(df_mantel, "mean_bray.16S", "fst")
  
  # Print results
  message("\n========== MANTEL TEST RESULTS ==========")
  message(paste("  16S Mantel r:", mantel_16S$r, "| p-value:", mantel_16S$p))
  message(paste("  ITS Mantel r:", mantel_ITS$r, "| p-value:", mantel_ITS$p))
  message(paste("  AMF Mantel r:", mantel_AMF$r, "| p-value:", mantel_AMF$p))
  
  
  # ADD RESULTS TO SUMMARY TABLE
  mantel_summary <- bind_rows(
    mantel_summary,
    tibble(
      taxa = TAXA,
      group = c("16S", "ITS", "AMF"),
      mantel.r = c(mantel_16S$r, mantel_ITS$r, mantel_AMF$r),
      mantel.p = c(mantel_16S$p, mantel_ITS$p, mantel_AMF$p)
    )
  )
  
  # Create plots for each dataset
  plots_16S <- create_fst_plots(df_mantel, "mean_bray.16S", "16S Bacteria", mantel_16S)
  plots_ITS <- create_fst_plots(df_mantel, "mean_bray.ITS", "ITS Fungi", mantel_ITS)
  plots_AMF <- create_fst_plots(df_mantel, "mean_bray.AMF", "ITS AMF", mantel_AMF)
  
  #### COMBINE PLOTS ####
  # Individual dataset plots (FST vs Bray + Ridge)
  combo_16S <- (plots_16S$with_self | plots_16S$no_self) / results_16S$ridge +
    plot_annotation(tag_levels = "A", title = paste0("16S Bacteria - ", TAXA)) +
    plot_layout(heights = c(1, 5))
  
  combo_ITS <- (plots_ITS$with_self | plots_ITS$no_self) / results_ITS$ridge +
    plot_annotation(tag_levels = "A", title = paste0("ITS Fungi - ", TAXA)) +
    plot_layout(heights = c(1, 5))
  
  combo_AMF <- (plots_AMF$with_self | plots_AMF$no_self) / results_AMF$ridge +
    plot_annotation(tag_levels = "A", title = paste0("ITS AMF - ", TAXA)) +
    plot_layout(heights = c(1, 5))
  
  # Combined FST vs Bray comparison (all three side-by-side)
  combo_fst_comparison <- (plots_16S$no_self | plots_ITS$no_self | plots_AMF$no_self) +
    plot_annotation(
      tag_levels = "A",
      title = paste0("FST vs Mean Bray-Curtis Distance Comparison - ", TAXA)
    )
  
  #### SAVE PLOTS ####
  message("\n========== SAVING PLOTS ==========")
  
  ggsave(paste0(
    OUTDIR_MANTEL_16S, "/0", counter, "-16S-fst-bray-no-self-", TAXA, ".jpg"), 
    plots_16S$no_self + 
      labs(title = paste0("16S - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_16S$r, ", p = ", mantel_16S$p)),
           width = 7, height = 6, dpi = 300)
  
  ggsave(paste0(
    OUTDIR_MANTEL_16S, "/0", counter, "-16S-fst-bray-no-self-", TAXA, ".svg"), 
    plots_16S$no_self + 
      labs(title = paste0("16S - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_16S$r, ", p = ", mantel_16S$p)),
    width = 7, height = 6, dpi = 300)
  
  ggsave(paste0(
    OUTDIR_MANTEL_ITS, "/0", counter, "-ITS-fst-bray-no-self-", TAXA, ".jpg"), 
    plots_ITS$no_self + 
      labs(title = paste0("ITS - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_ITS$r, ", p = ", mantel_ITS$p)),
    width = 7, height = 6, dpi = 300)
  ggsave(paste0(
    OUTDIR_MANTEL_ITS, "/0", counter, "-ITS-fst-bray-no-self-", TAXA, ".svg"), 
    plots_16S$no_self + 
      labs(title = paste0("16S - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_ITS$r, ", p = ", mantel_ITS$p)),
    width = 7, height = 6, dpi = 300)
  
  ggsave(paste0(
    OUTDIR_MANTEL_AMF, "/0", counter, "-AMF-fst-bray-no-self-", TAXA, ".jpg"), 
    plots_AMF$no_self + 
      labs(title = paste0("AMF - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_AMF$r, ", p = ", mantel_AMF$p)),
    width = 7, height = 6, dpi = 300)
  ggsave(paste0(
    OUTDIR_MANTEL_AMF, "/0", counter, "-AMF-fst-bray-no-self-", TAXA, ".svg"), 
    plots_AMF$no_self + 
      labs(title = paste0("AMF - ",TAXA), 
           subtitle = paste0("Mantel r = ", mantel_AMF$r, ", p = ", mantel_AMF$p)),
    width = 7, height = 6, dpi = 300)
  
  ggsave(paste0(OUTDIR_RIDGE_16S, "/0", counter, "-16S-fst-bray-ridge-", TAXA, 
                ".jpg"), combo_16S, width = 10, height = 12, dpi = 300)
  ggsave(paste0(OUTDIR_RIDGE_16S, "/0", counter, "-16S-fst-bray-ridge-", TAXA, 
                ".svg"), combo_16S, width = 10, height = 12, dpi = 300)
  
  ggsave(paste0(OUTDIR_RIDGE_ITS, "/0", counter, "-ITS-fst-bray-ridge-", TAXA, 
                ".jpg"), combo_ITS, width = 10, height = 12, dpi = 300)
  ggsave(paste0(OUTDIR_RIDGE_ITS, "/0", counter, "-ITS-fst-bray-ridge-", TAXA, 
                ".svg"), combo_ITS, width = 10, height = 12, dpi = 300)
  
  ggsave(paste0(OUTDIR_RIDGE_AMF, "/0", counter, "-AMF-fst-bray-ridge-", TAXA, 
                ".jpg"), combo_AMF, width = 10, height = 12, dpi = 300)
  ggsave(paste0(OUTDIR_RIDGE_AMF, "/0", counter, "-AMF-fst-bray-ridge-", TAXA, 
                ".svg"), combo_AMF, width = 10, height = 12, dpi = 300)
  
  ggsave(paste0(OUTDIR_MANTEL_TOGETHER, "/0", counter, "-mantel-comparison-", 
                TAXA, ".jpg"), combo_fst_comparison, width = 15, height = 5, 
         dpi = 300)
  ggsave(paste0(OUTDIR_MANTEL_TOGETHER, "/0", counter, "-mantel-comparison-", 
                TAXA, ".svg"), combo_fst_comparison, width = 15, height = 5, 
         dpi = 300)
  
  
  message("\n========== ALL ANALYSES COMPLETE ==========")

}
  
#AFTER the loop, save the summary table
message("\n========== SAVING MANTEL SUMMARY TABLE ==========")
print(mantel_summary, n = Inf)


mantel_wide <- mantel_summary %>%
  pivot_wider(
    names_from = group,           # Column to spread into new columns
    values_from = c(mantel.r, mantel.p)  # Values to fill those columns
  ) %>% 
  select(taxa, mantel.r_16S, mantel.p_16S, mantel.r_ITS, mantel.p_ITS, 
         mantel.r_AMF, mantel.p_AMF)

write_xlsx(mantel_wide, "plots/mantel-overview/mantel_summary.xlsx")
