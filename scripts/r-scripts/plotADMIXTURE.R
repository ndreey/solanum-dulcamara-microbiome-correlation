#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
admix_prefix <- args[1]
samples <- args[2]
out_dir <- args[3]
K <- as.numeric(args[4])

# Create output directory if it doesn't exist
output_dir <- out_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===== Arguments =====
POP_ORDER <- c("Kämpinge","Malmo","Alnarp_pond", "Lomma_beach", "Habo_Ljung", "Lund", "Borgeby", "Tvedöra")
Q_FILE   <- sprintf(paste0(admix_prefix,".%d.Q"), K)
FAM_FILE <- samples
META_FILE <- "doc/S_dulcamara_popmap_1.txt"       # has ID, Population, etc.

# ---- 1) read Q matrix ----
df_q <- read_delim(Q_FILE, delim = " ", col_names = FALSE)
colnames(df_q) <- paste0("Q", seq_len(K))

# ---- 2) read IDs + metadata ----
ids  <- read_delim(FAM_FILE,  delim = "\t", col_names = "ID")
meta <- read_delim(META_FILE, delim = "\t", col_names = TRUE)

# ---- 3) attach IDs + metadata (wide form for ordering) ----
wide <- df_q %>%
  mutate(ID = ids$ID) %>%
  left_join(meta, by = "ID") %>%
  mutate(Population = factor(Population, levels = POP_ORDER))

# ---- 4) order individuals within each population by Q2 (desc) ----
SORT_CLUSTER <- "Q2"   # change to "Q1" etc. if you prefer
ordered_ids <- wide %>%
  arrange(Population) %>%
  group_by(Population) %>%
  arrange(desc(.data[[SORT_CLUSTER]]), .by_group = TRUE) %>%
  ungroup() %>%
  pull(ID) %>%
  unique()

# ---- 5) long format for plotting ----
long_df <- wide %>%
  pivot_longer(starts_with("Q"), names_to = "Cluster", values_to = "Proportion") %>%
  mutate(
    ID = factor(ID, levels = rev(ordered_ids)),
    Population = factor(Population, levels = POP_ORDER)
  )

# ---- 6) main plot (no grey facet strip) ----
p_admix <- ggplot(long_df, aes(ID, Proportion, fill = Cluster)) +
  geom_col(width = 0.99) +
  facet_grid(~ Population, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),  # removes grey rectangle
    strip.text = element_text(size = 14)
  )

label_plot <- ggdraw() + draw_label(paste0("K = ", K), angle = 270, size = 14)
final_plot <- plot_grid(p_admix, label_plot, ncol = 2, rel_widths = c(1, 0.03))

ggsave(file.path(output_dir, paste0("admixture-",K,".png")), final_plot, 
    width = 14, height = 4, dpi = 300)