#!/usr/bin/env Rscript

library(tidyverse)
library(scales)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
frq_file <- args[1]
idepth_file <- args[2]
imiss_file <- args[3]
lmiss_file <- args[4]
out_dir <- args[5]

# Create output directory if it doesn't exist
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Loading files")

imiss <- read_delim(imiss_file, delim = "\t", col_names = T)
idepth <- read_delim(idepth_file, delim = "\t", col_names = T)
frq <- read_delim(frq_file, delim = "\t", skip = 1, 
                  col_names = c("CHROM", "POS", "N_A", 
                                "N_C", "FreqA", "FreqB"))

lmiss <- read_delim(lmiss_file, delim = "\t", col_names = T)

# find minor allele frequency
frq$maf <- frq %>% select(FreqA, FreqB) %>% apply(1, function(z) min(z))


meta_data <- read_delim("doc/S_dulcamara_popmap_1.txt", delim = "\t", col_names = T)

# Join metadata
df_imiss <- imiss %>%
  inner_join(meta_data, by = c("INDV" = "ID"))

df_idepth <- idepth %>%
  inner_join(meta_data, by = c("INDV" = "ID"))

# Boxplot of F_MISS per individual
p_imiss <- ggplot(df_imiss, aes(x = Population, y = F_MISS)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_hline(yintercept=0.3, linetype="dashed", color = "gray") +
labs(x = "Population", y = "F_MISS", title = "Individual F_MISS by Population") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file.path(output_dir, "imiss_boxplot.png"), p_imiss, width = 8, height = 6)

# Boxplot of MEAN_DEPTH
p_idepth <- ggplot(df_idepth, aes(x = Population, y = MEAN_DEPTH)) +
  geom_boxplot(outliers=F) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Population", y = "MEAN_DEPTH", title = "Mean Depth by Population") +
  theme_bw() +
  theme(panel.grid = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file.path(output_dir, "idepth_boxplot.png"), p_idepth, width = 8, height = 6)

# Density plot of locus missingness
p_lmiss <- ggplot(lmiss, aes(x = F_MISS)) +
  geom_density(fill = "red", alpha = 0.5) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  labs(x = "Locus Missingness (F_MISS)", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file.path(output_dir, "lmiss_density.png"), p_lmiss, width = 8, height = 6)

# Density plot of MAF
p_maf <- ggplot(frq, aes(x = maf)) +
  geom_density(fill = "blue", alpha = 0.5) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  labs(title = "Density Plot of Minor Allele Frequency (MAF)",
       x = "Minor Allele Frequency (MAF)", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file.path(output_dir, "maf_density.png"), p_maf, width = 8, height = 6)


# Define thresholds ---
thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05)

# Clean MAF vector (drop NAs) ---
maf <- frq$maf
maf <- maf[!is.na(maf)]

# Total SNPs with a computable MAF ---
total <- length(maf)

# Count how many SNPs would be removed at each threshold
n_removed <- numeric(length(thresholds))
for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  n_removed[i] <- sum(maf <= t)
}

# Build results data frame ---
lost <- data.frame(
  threshold    = thresholds,
  n_removed    = n_removed,
  n_kept       = total - n_removed,
  total        = total,
  prop_removed = n_removed / total
)

p_lost <- ggplot(lost, aes(x = threshold, y = n_removed)) +
  geom_line() +
  geom_point(size = 2) +
  geom_text(aes(label = n_removed), vjust = -0.6, size = 3) +
  scale_x_continuous(breaks = thresholds) +
  labs(
    title = "SNPs removed vs. MAF threshold",
    x = "MAF threshold (keeping MAF > threshold)",
    y = "Number of SNPs removed"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file.path(output_dir, "maf_threshold_scree.png"), p_lost, width = 8, height = 6)