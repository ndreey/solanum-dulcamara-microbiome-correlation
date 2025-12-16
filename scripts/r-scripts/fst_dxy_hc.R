library(tidyverse)
library(scales)
library(ggpubr)
library(patchwork)
library(paletteer)
library(pheatmap)

set.seed(1337)

# Load in data
df_piawka <- read_delim("04-piawka/GOR_solanum.piawka-30miss.hudson.tsv", delim = "\t",
                        col_names = T)
pop_map <- read_delim("doc/uniq_pop_meta.tsv", delim = "\t",
                      col_names = T)

# Wrangle the dataframes to get pair names, fst and geodiffs.
df_fst <- df_piawka %>% 
  select(pop1, pop2, metric, value) %>% 
  filter(!metric %in% c("pi", "Dxy")) %>% 
  rename(fst = value) %>% 
  select(pop1, pop2, fst) %>% 
  # Join for pop1
  left_join(pop_map %>% select(pop, pop.id), by = c("pop1" = "pop")) %>%
  rename(pop1.id = pop.id) %>%
  # Join for pop2  
  left_join(pop_map %>%  select(pop, pop.id), by = c("pop2" = "pop")) %>%
  rename(pop2.id = pop.id) %>% 
  #mutate(pair = paste0(pop1.id,"-",pop2.id)) %>% 
  select(pop1.id, pop2.id, fst, -pop1, -pop2)

# Wrangle the dataframes to get pair names, fst and geodiffs.
df_dxy <- df_piawka %>% 
  select(pop1, pop2, metric, value) %>% 
  filter(!metric %in% c("pi", "Fst_HUD")) %>% 
  rename(dxy = value) %>% 
  select(pop1, pop2, dxy) %>% 
  # Join for pop1
  left_join(pop_map %>% select(pop, pop.id), by = c("pop1" = "pop")) %>%
  rename(pop1.id = pop.id) %>%
  # Join for pop2  
  left_join(pop_map %>%  select(pop, pop.id), by = c("pop2" = "pop")) %>%
  rename(pop2.id = pop.id) %>% 
  #mutate(pair = paste0(pop1.id,"-",pop2.id)) %>% 
  select(pop1.id, pop2.id, dxy, -pop1, -pop2)

# Make symmetric matrix
df_sym <- df_fst %>% 
  bind_rows(
    df_fst,
    df_fst %>% rename(pop1.id = pop2.id, pop2.id = pop1.id)
  ) %>% 
  bind_rows(
    tibble(
      pop1.id = unique(c(df_fst$pop1.id, df_fst$pop2.id)),
      pop2.id = unique(c(df_fst$pop1.id, df_fst$pop2.id)),
      fst = 0
    )
  ) %>% 
  distinct() %>% 
  arrange(pop1.id, pop2.id)

# Make symmetric matrix
df_sym_dxy <- df_dxy %>% 
  bind_rows(
    df_dxy,
    df_dxy %>% rename(pop1.id = pop2.id, pop2.id = pop1.id)
  ) %>% 
  bind_rows(
    tibble(
      pop1.id = unique(c(df_dxy$pop1.id, df_dxy$pop2.id)),
      pop2.id = unique(c(df_dxy$pop1.id, df_dxy$pop2.id)),
      dxy = 0
    )
  ) %>% 
  distinct() %>% 
  arrange(pop1.id, pop2.id)

######## Plot fst

# Convert to matrix
fst_matrix <- df_sym %>% 
  pivot_wider(names_from = pop2.id, values_from = fst) %>% 
  column_to_rownames("pop1.id") %>% 
  as.matrix()


# Create annotation dataframe from your pop_map
# The rownames MUST match your matrix row/column names
annotation_df <- pop_map %>% 
  select(pop.id, habitat, humidity) %>% 
  column_to_rownames("pop.id")

# Check that rownames match
all(rownames(annotation_df) %in% rownames(fst_matrix))  # Should be TRUE

# Define colors for each annotation
ann_colors <- list(
  habitat = c(
    forest = "chocolate3",
    urban = "maroon3", 
    rural = "mediumpurple",
    beach = "aquamarine4"
  ),
  humidity = c(
    wet = "slateblue",
    dry = "tan"
  )
)


pheatmap(
  fst_matrix,
  cellwidth=20,
  cellheight=20,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  treeheight_col = 0,
  
  # Add annotations (same for rows and columns since it's symmetric)
  annotation_row = annotation_df,
  #annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_legend = T,
  
  fontsize = 12,
  legend = TRUE, border_color = NA,
  color = rev(paletteer_c("grDevices::Blues", n = 100)),
  display_numbers = F,
  angle_col = 45,
  filename = "plots/fst_small_heatmap.png")



# Convert to matrix
dxy_matrix <- df_sym_dxy %>% 
  pivot_wider(names_from = pop2.id, values_from = dxy) %>% 
  column_to_rownames("pop1.id") %>% 
  as.matrix()

pheatmap(
  dxy_matrix,
  cellwidth=20,
  cellheight=20,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  treeheight_col = 0,
  
  # Add annotations (same for rows and columns since it's symmetric)
  annotation_row = annotation_df,
  #annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_legend = T,
  
  fontsize = 12,
  legend = TRUE, border_color = NA,
  color = rev(paletteer_c("grDevices::Reds", n = 100)),
  display_numbers = F,
  angle_col = 45,
  filename = "plots/dxy_small_heatmap.png")











# Create composite matrix
# Lower triangle = FST, Upper triangle = Dxy, Diagonal = 0

composite_matrix <- fst_matrix  # Start with FST matrix

# Replace upper triangle with Dxy values
composite_matrix[upper.tri(composite_matrix)] <- dxy_matrix[upper.tri(dxy_matrix)]

# Visual check
print(composite_matrix, digits = 3)


#Convert matrix to data frame with rownames as a column
composite_df <- composite_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("Population")


library(writexl)
# Write to Excel
write_xlsx(composite_df, "composite_fst_dxy_matrix.xlsx")










