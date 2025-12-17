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
# This script will create "master mantel data-frames" to build upon
# to make more mantel tests. Easiest, is to use the pop_map to map the
# pairwise comparison measures.
#

source("scripts/r-scripts/master-data-frame-functions.R")

set.seed(1337)

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

# Load in km differences
geo <- read_delim("doc/geo_distance.tsv", delim = "\t", col_names = T) %>% 
  separate(pair, into = c("pop1.id", "pop2.id"), sep = "-", remove = FALSE)

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

# Make geo symmetric
geo_sym <- geo %>%
  bind_rows(
    geo %>% rename(pop1.id = pop2.id, pop2.id = pop1.id)
  # Add self-comparisons with 0 distance
  ) %>% 
  bind_rows(
    tibble(
      pop1.id = unique(c(geo$pop1.id, geo$pop2.id)),
      pop2.id = unique(c(geo$pop1.id, geo$pop2.id)),
      km.diff = 0
    )
  ) %>%
  mutate(pair = paste0(pop1.id, "-", pop2.id))

#### Loop through taxonomic levels #####
# List to store master frames per TAXA
master_list <- list()

for (TAXA in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  message(paste("Processing master dataframe for", TAXA))

  # Load abundance data
  abund_16S <- read_delim(paste0("data/16S_abundance/", TAXA, "_abund.csv"),
                          delim = ",", col_names = T) %>% 
    rename("taxon" = "...1") %>% 
    filter(taxon != "k__Archaea") %>% 
    # Filter to only retain S. dulcamara samples
    column_to_rownames("taxon") %>%
    select(all_of(meta_dulcamara$sampleID))
  
  abund_ITS <- read_delim(paste0("data/ITS_abundance/", TAXA, "_abund.csv"),
                          delim = ",", col_names = T) %>% 
    rename("taxon" = "...1") %>% 
    # Filter to only retain S. dulcamara samples
    column_to_rownames("taxon") %>%
    select(all_of(meta_dulcamara$sampleID))
  
  
  # AMF data is in different format, wrangle it to match
  abund_AMF <- read_delim(paste0("data/AMF_abundance/", TAXA, "_abund.tsv"), 
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
  
  
  results_16S <- get_bray_curtis(abund_16S, "16S", meta_dulcamara)
  results_ITS <- get_bray_curtis(abund_ITS, "ITS", meta_dulcamara)
  results_AMF <- get_bray_curtis(abund_AMF, "AMF", meta_AMF)
  
  #### MERGE DATA ####
  df_mantel <- df_sym %>% 
    inner_join(results_AMF, by = "pair") %>% 
    rename(mean_bray.AMF = mean_bray) %>% 
    inner_join(results_ITS, by = "pair") %>% 
    rename(mean_bray.ITS = mean_bray) %>% 
    inner_join(results_16S, by = "pair") %>% 
    rename(mean_bray.16S = mean_bray) %>% 
    select(pair, pop1.id, pop2.id, fst, starts_with("mean_bray")) %>% 
    inner_join(geo_sym %>% select(pair, km.diff), by = "pair") %>% 
    arrange(pop1.id, pop2.id)
  
  # Store in list
  master_list[[TAXA]] <- df_mantel
  message(paste0("Saved mantel dataframe of ", TAXA))
  
}

# All mantel dataframes stored in list at their taxonomic level.
master_frames <- master_list

# Save the object to an .Rds file 
saveRDS(master_frames, file = "data/master_frames.Rds") 

# Combine all taxonomic levels into one dataframe
df_all <- bind_rows(
  master_frames$Phylum %>% mutate(taxon_level = "Phylum"),
  master_frames$Class %>% mutate(taxon_level = "Class"),
  master_frames$Order %>% mutate(taxon_level = "Order"),
  master_frames$Family %>% mutate(taxon_level = "Family"),
  master_frames$Genus %>% mutate(taxon_level = "Genus"),
  master_frames$Species %>% mutate(taxon_level = "Species")
)

# Save as .csv and .Rds
write_csv(df_all, "data/mantel_all_taxa.csv")
saveRDS(df_all, file = "data/mantel_all_taxa.Rds")

# Later, you can load it back with: 
#master_frames <- readRDS("data/master_frames.Rds")
#df_all <- readRDS("data/mantel_all_taxa.Rds")


#### Exampel usage ####
# They can be separated out like this.
df_phylum <- master_frames$Phylum
df_class <- master_frames$Class
df_order <- master_frames$Order
df_family <- master_frames$Family
df_genus <- master_frames$Genus
df_species <- master_frames$Species

##### Save the Master Dataframes as .csv ######
write_csv(df_phylum, "data/mantel_phylum.csv")
write_csv(df_class, "data/mantel_class.csv")
write_csv(df_order, "data/mantel_order.csv")
write_csv(df_family, "data/mantel_family.csv")
write_csv(df_genus, "data/mantel_genus.csv")
write_csv(df_species, "data/mantel_species.csv")




# Prepare data for plotting (exclude duplicated pairs and self-comparisons)
df_plot <- df_family %>%
  filter(pop1.id < pop2.id)

# Simple scatter plot
ggplot(df_plot, aes(x = fst, y = km.diff)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = lm, se = FALSE, color = "red", linetype = "dashed") +
  theme_pubr()

# Run Mantel tests
mantel_fam <- run_mantel(df_family, "mean_bray.16S", "fst")
cat("Mantel r =", mantel_fam$r, "\n")
cat("p-value =", mantel_fam$p, "\n")
  
# Summarize across beta-div. 
beta_div <- c("mean_bray.16S", "mean_bray.ITS", "mean_bray.AMF")
comparison_results <- data.frame()

# Re-run this for loop a couple times and see how correlation is consistent
# but p-value varies (even though we have set seed).
for (bd in beta_div) {
  result <- run_mantel(df_family, bd, "fst")
  
  comparison_results <- rbind(
    comparison_results,
    data.frame(
      microbiome = bd,
      mantel_r = result$r,
      p_value = result$p,
      significant = result$p < 0.05
    )
  )
}

cat("\n=== Comparing Different Microbiome Types ===\n")
print(comparison_results)
#> print(comparison_results)
#microbiome mantel_r p_value significant
#1  mean_bray.16S   0.5100   0.051       FALSE
#2  mean_bray.ITS   0.2115   0.232       FALSE
#3  mean_bray.AMF   0.0874   0.399       FALSE
#4  mean_bray.16S   0.5100   0.039        TRUE
#5  mean_bray.ITS   0.2115   0.234       FALSE
#6  mean_bray.AMF   0.0874   0.385       FALSE
#7  mean_bray.16S   0.5100   0.045        TRUE
#8  mean_bray.ITS   0.2115   0.233       FALSE
#9  mean_bray.AMF   0.0874   0.367       FALSE
#10 mean_bray.16S   0.5100   0.046        TRUE
#11 mean_bray.ITS   0.2115   0.244       FALSE
#12 mean_bray.AMF   0.0874   0.390       FALSE
#13 mean_bray.16S   0.5100   0.040        TRUE
#14 mean_bray.ITS   0.2115   0.257       FALSE
#15 mean_bray.AMF   0.0874   0.418       FALSE
#16 mean_bray.16S   0.5100   0.050       FALSE
#17 mean_bray.ITS   0.2115   0.253       FALSE
#18 mean_bray.AMF   0.0874   0.393       FALSE



##### Save the Master Dataframes as .csv ######
write_csv(df_phylum, "data/mantel_phylum.csv")
write_csv(df_class, "data/mantel_class.csv")
write_csv(df_order, "data/mantel_order.csv")
write_csv(df_family, "data/mantel_family.csv")
write_csv(df_genus, "data/mantel_genus.csv")
write_csv(df_species, "data/mantel_species.csv")