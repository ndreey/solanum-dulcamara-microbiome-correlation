library(tidyverse)
library(dada2)
library(RcppParallel)
library(phyloseq)
library(readxl)

message("========== LOADING DADA2 OBJECTS ==========")
# Load all RDS objects
out <- readRDS("amplicon_rds/ITS/out.rds")
seqtab <- readRDS("amplicon_rds/ITS/seqtab.rds")
seqtab.nochim <- readRDS("amplicon_rds/ITS/seqtab.nochim.rds")
taxa <- readRDS("amplicon_rds/ITS/taxa.rds")

message(paste("Samples:", nrow(seqtab.nochim), "| ASVs:", ncol(seqtab.nochim)))
message(paste("Reads retained:", round(sum(seqtab.nochim)/sum(seqtab), 3)))

message("\n========== CREATING OTU TABLE WITH TAXONOMY ==========")

# Load and prepare metadata
metadata <- read_excel("doc/raw-meta.xlsx") %>% as.data.frame()
rownames(metadata) <- metadata$sampleID

# Load common data once (outside loop)
meta_dulcamara <- read_delim("doc/meta-microeco.tsv", delim = "\t", col_names = T) %>% 
  filter(spp == "Solanum_dulcamara")

# Fix sample names (ITS → D)
rownames(seqtab.nochim) <- str_replace(rownames(seqtab.nochim), "ITS", "D")

# Create phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
  sample_data(metadata), 
  tax_table(taxa)
)

# Extract and merge tables
otu_df_raw <- otu_table(ps) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  left_join(
    tax_table(ps) %>%
      as.data.frame() %>%
      rownames_to_column("ASV_ID") %>%
      mutate(
        Kingdom = ifelse(is.na(Kingdom), "k__", Kingdom),
        Phylum = ifelse(is.na(Phylum), "p__", Phylum),
        Class = ifelse(is.na(Class), "c__", Class),
        Order = ifelse(is.na(Order), "o__", Order),
        Family = ifelse(is.na(Family), "f__", Family),
        Genus = ifelse(is.na(Genus), "g__", Genus),
        Species = ifelse(is.na(Species), "s__", Species)
      ),
    by = "ASV_ID"
  )

# ========== SAVE RAW UNFILTERED OTU TABLE ==========
message("\n========== SAVING RAW OTU TABLE ==========")

# Create taxonomy string
otu_raw_export <- otu_df_raw %>%
  mutate(
    taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "|"),
    .after = "ASV_ID"
  ) %>%
  select(taxonomy, starts_with("D"))

# Save
write_tsv(otu_raw_export, "otu_tables/ITS-raw-all-ASVs.tsv")
message(paste("Saved raw OTU table:", nrow(otu_raw_export), "ASVs"))

# ========== FILTER TO DULCAMARA AND NO ZERO TAXA ==========
message("\n========== FILTER TO DULCAMARA AND NO ZERO TAXA ==========")

message("Before filtering:")
message(paste("  Rows:", nrow(otu_df_raw)))
message(paste("  Columns:", ncol(otu_df_raw)))

otu_df <- otu_df_raw %>%
  select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
         all_of(meta_dulcamara$sampleID)) %>%
  # Keep only rows where sum of D columns > 0
  filter(rowSums(select(., starts_with("D"))) > 0)

message("\nAfter filtering:")
message(paste("  Rows:", nrow(otu_df)))
message(paste("  Columns:", ncol(otu_df)))
message(paste("  Samples kept:", length(meta_dulcamara$sampleID)))
message(paste("  ASVs removed:", nrow(otu_df_raw) - nrow(otu_df)))



message("\n========== CREATING TAXONOMIC TABLES ==========")
# Generic function
create_tax_tables <- function(otu_df, filter_expr, prefix, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  otu_filtered <- otu_df %>% filter({{filter_expr}})
  message(paste("Filtered to", nrow(otu_filtered), prefix, "ASVs"))
  
  tax_levels <- list(
    phylum = c("Kingdom", "Phylum"),
    class = c("Kingdom", "Phylum", "Class"),
    order = c("Kingdom", "Phylum", "Class", "Order"),
    family = c("Kingdom", "Phylum", "Class", "Order", "Family"),
    genus = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
    species = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
  
  for (level in names(tax_levels)) {
    levels_to_keep <- tax_levels[[level]]
    
    table <- otu_filtered %>%
      group_by(across(all_of(levels_to_keep))) %>%
      summarise(across(starts_with("D"), sum), .groups = "drop") %>%
      mutate(taxonomy = paste(!!!syms(levels_to_keep), sep = "|"), .before = 1) %>%
      select(-all_of(levels_to_keep))
    
    filename <- file.path(output_dir, paste0(prefix, "-", level, ".tsv"))
    write_tsv(table, filename)
    message(paste("  ", level, ":", nrow(table), "taxa →", basename(filename)))
  }
  
  message("All tables created!\n")
}

# Create Fungi tables
create_tax_tables(otu_df, Kingdom == "k__Fungi", "ITS-fungi", "otu_tables/ITS-fungi")

# Create AMF tables
create_tax_tables(otu_df, Phylum == "p__Glomeromycota", "ITS-AMF", "otu_tables/ITS-AMF")
