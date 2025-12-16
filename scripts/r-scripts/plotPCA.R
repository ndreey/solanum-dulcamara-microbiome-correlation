#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(ggrepel)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
eigenval_path <- args[1]
eigenvec_path <- args[2]
out_dir <- args[3]

# Create output directory if it doesn't exist
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===== LOAD META DATA =====
meta_data <- read_delim("doc/S_dulcamara_popmap_1.txt", delim = "\t", col_names = T)

# ===== LOAD PCA DATA =====
message("Loading files")
eigenval <- scan(eigenval_path)
pca <- read_delim(eigenvec_path, col_names = FALSE)[, -1]
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))
percent_var <- eigenval / sum(eigenval) * 100
pve <- data.frame(PC = 1:length(percent_var), pve = percent_var)

pca <- pca %>%
  mutate(across(matches("^PC\\d+$"), ~ as.numeric(.)))
  
# ===== Merge Meta Data =====
message("Merge with meta data")
pca <- pca %>%
  inner_join(meta_data, by = c("ID" = "ID"))

# make sure these are discrete
pca <- pca %>% mutate(Population = factor(Population),
                      Habitat        = factor(Habitat))
# ===== PLOTS =====
## Scree Plot
message("Plotting screeplot")
p_scree <- ggplot(pve, aes(x = PC, y = pve)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f%%", pve)), vjust = -0.5, size = 2.5) +
  scale_x_continuous(breaks = pve$PC, limits = c(0.5, max(pve$PC) + 0.5)) +
  labs(title = "Scree Plot of PCA",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(file.path(output_dir, "PCA_scree.png"), p_scree, width = 10, height = 6)



pc12 <- pca %>%
  ggplot(aes(PC1, PC2, fill = Population, shape = Habitat)) +
  geom_point(size = 4, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  scale_fill_brewer(palette = "Paired") +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black")),
    shape = guide_legend(override.aes = list(fill  = "grey80")) # optional
  ) +
  geom_text_repel(
    aes(label = ID),
    size = 3,
    max.overlaps = 5,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 123
    # if you want labels colored by population, add: colour = Population
  ) +
  scale_x_continuous(name = paste0("PC1 (", round(percent_var[1], 2), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(percent_var[2], 2), "%)")) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", size = 1.2))
ggsave(file.path(output_dir, "PCA_1vs2.png"), pc12, width = 10, height = 6)

# PC3 vs PC4 (filled shapes)
pc34 <- pca %>%
  ggplot(aes(PC3, PC4, fill = Population, shape = Habitat)) +
  geom_point(size = 4, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_brewer(palette = "Paired") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black")),
    shape = guide_legend(override.aes = list(fill  = "grey80")) # optional
  ) +
    geom_text_repel(
    aes(label = ID),
    size = 3,
    max.overlaps = 5,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 123
  ) +
  scale_x_continuous(name = paste0("PC3 (", round(percent_var[3], 2), "%)")) +
  scale_y_continuous(name = paste0("PC4 (", round(percent_var[4], 2), "%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "right",
    panel.border = element_rect(color = "black", size = 1.2)
  )
ggsave(file.path(output_dir, "PCA_3vs4.png"), pc34, width = 10, height = 6)

# PC3 vs PC4 (filled shapes)
pc56 <- pca %>%
  ggplot(aes(PC5, PC6, fill = Population, shape = Habitat)) +
  geom_point(size = 4, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_brewer(palette = "Paired") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", alpha = 0.5) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, colour = "black")),
    shape = guide_legend(override.aes = list(fill  = "grey80")) # optional
  ) +
    geom_text_repel(
    aes(label = ID),
    size = 3,
    max.overlaps = 5,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 123
  ) +
  scale_x_continuous(name = paste0("PC3 (", round(percent_var[5], 2), "%)")) +
  scale_y_continuous(name = paste0("PC4 (", round(percent_var[6], 2), "%)")) +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "right",
    panel.border = element_rect(color = "black", size = 1.2)
  )
ggsave(file.path(output_dir, "PCA_5vs6.png"), pc56, width = 10, height = 6)


message("Rscript done")