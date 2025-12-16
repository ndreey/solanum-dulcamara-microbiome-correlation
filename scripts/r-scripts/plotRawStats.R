#!/usr/bin/env Rscript

library(tidyverse)
library(scales)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
stats_tsv <- args[1]
out_dir <- args[2]

# Create output directory if it doesn't exist
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Loading file: ", stats_tsv)
df <- read_delim(stats_tsv, delim = "\t", col_names = T)


message("Plotting AN")
# Plot density plot  for AN
p_AN <- df %>%
  ggplot(aes(x = AN)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(breaks = pretty_breaks(n=15), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AN_density.png"),
       p_AN, width = 8, height = 6)

message("Plotting AC")
# Plot density plot  for AC
p_AC <- df %>%
  ggplot(aes(x = AC)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(breaks = pretty_breaks(n=15), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AC_density.png"),
       p_AC, width = 8, height = 6)

message("Plotting AF")
# Plot density plot  for AF
p_AF <- df %>%
  ggplot(aes(x = AF)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(breaks = pretty_breaks(n=15), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "AF_density.png"),
       p_AF, width = 8, height = 6)

message("Plotting ExcHet")
# Plot density plot  for ExcHet
p_ExcHet <- df %>%
  ggplot(aes(x = ExcHet)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(breaks = pretty_breaks(n=15), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "ExcHet_density.png"),
       p_ExcHet, width = 8, height = 6)



message("Plotting F_MISSING")
# Plot density plot  for ReadPosRankSum
p_FMISS <- df %>% 
  ggplot(aes(x = F_MISSING)) +
  geom_density(alpha = 0.5, fill = "darkolivegreen") +
  scale_x_continuous(breaks = pretty_breaks(n=15), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank())

ggsave(file.path(output_dir, "F_MISSING.png"),
       p_FMISS, width = 8, height = 6)

message("Plotting complete")