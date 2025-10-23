library(circlize)
library(ComplexHeatmap)

# derep
greening_hits_best <- greening_hits_filt %>%
  group_by(sample, gene_id) %>%
  slice_max(pident, n = 1) %>%
  ungroup()

heatmap_data <- greening_hits_best %>%
  group_by(sample, database_clean, category, original_n_reads) %>%
  dplyr::summarise(count = n(), length = mean(length), .groups = "drop") %>%
  mutate(rpm = count * 1e9/original_n_reads/length) %>%
  select(-count, -original_n_reads, -length) %>%
  pivot_wider(names_from = sample, values_from = rpm, values_fill = 0)

# create matrix
mat <- as.matrix(heatmap_data[ , c("A-K", "B-K", "E-K", "G-K", "H-K")])
rownames(mat) <- heatmap_data$database_clean

# order rows within each category
heatmap_data <- heatmap_data %>%
  group_by(category) %>%
  dplyr::arrange(desc(rowSums(dplyr::across(all_of(c("A-K","B-K","E-K","G-K","H-K")))))) %>%
  ungroup()
mat <- mat[heatmap_data$database_clean, ]

# color categories
categories <- unique(heatmap_data$category)

category_colors <- c("Aerobic respiration" = "#D81B60",
                     "Alternative e acceptor"= "#1E88E5",
                     "Alternative e donor"= "#10819D",
                     "Carbon fixation"  = "#DC2704",
                     "Nitrogen metabolism" = "#004D40",
                     "Trace gas metabolism" = "#766B4B",
                     "Phototrophy" = "#4DBE64",
                     "Sulfur metabolism" = "#FFC107",
                     "Other" = "grey85"
)

row_ha <- rowAnnotation(
  Category = heatmap_data$category,
  col = list(Category = category_colors),
  show_annotation_name = FALSE
)

# plot heatmap
genes_heatmap_RPM <- Heatmap(
  mat,
  name = "RPKM",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_split = heatmap_data$category,
  row_title = NULL,
  left_annotation = row_ha,
  col = colorRamp2(c(0, max(mat)), c("grey95", "#453781ff")),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = FALSE,
  row_names_max_width = unit(15, "cm"),
  column_names_max_height = unit(5, "cm"),
  column_names_centered = TRUE
  #  cell_fun = function(j, i, x, y, width, height, fill) { # this is to show actual values on the map
  # grid.text(
  #   round(mat[i, j], 1),
  #   x = x, y = y,
  #   gp = gpar(fontsize = 6, col = "black")
  # )
  # }
)

# save as .pdf
pdf("~.pdf", width = 7, height = 10) 
draw(genes_heatmap_RPM)
dev.off()