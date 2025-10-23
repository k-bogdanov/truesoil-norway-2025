# Prepare data
hits_and_properties <- left_join(greening_hits_best,
                                 Norway_metadata %>%
                                   select(-SampleID, -Treatment) %>%
                                   filter(!str_detect(PlotID, "R")) %>%
                                   dplyr::rename("sample" = PlotID),
                                 by = "sample"
) %>%
  group_by(sample, database_clean, category, original_n_reads, Distance_from_top_m, Bulk_C_mg_g_dw, Bulk_N_mg_g_dw, pH_H2O, Clay, Silt, Sand,
           MAOC_mg_g_dw, POC_mg_g_dw, Bulk_density_g_cm3, Particle_density_g_cm3, Volumetric_soil_water_content,
           Average_litter_derived_CO2_mg_g_dw, Average_soil_derived_CO2_mg_g_dw,
           Average_total_CO2_mg_g_dw) %>%
  dplyr::summarise(count = n(), length = mean(length), .groups = "drop") %>%
  mutate(rpm = count * 1e9/original_n_reads/length) %>%
  dplyr::rename("SOC (mg/g dw)" = Bulk_C_mg_g_dw,
                "N (mg/g dw)" = Bulk_N_mg_g_dw,
                "Distance from top (m)" = Distance_from_top_m,
                "pH" = pH_H2O,
                "Clay (%)" = Clay,
                "Silt (%)" = Silt,
                "Sand (%)" = Sand,
                "Vol. water content (%)" = Volumetric_soil_water_content,
                "MAOC (mg/g dw)" = MAOC_mg_g_dw,
                "POC (mg/g dw)" = POC_mg_g_dw,
                "Bulk density (g/cm^3)" = Bulk_density_g_cm3,
                "Particle density (g/cm^3)" = Particle_density_g_cm3,
                "Litter-derived CO2 mg/g dw" = Average_litter_derived_CO2_mg_g_dw,
                "Soil-derived CO2 mg/g dw" = Average_soil_derived_CO2_mg_g_dw,
                "Total CO2 mg/g dw" = Average_total_CO2_mg_g_dw)

soil_variable <- hits_and_properties %>%
  select(-sample, -database_clean, -category, -original_n_reads, -count, -length, -rpm)

soil_variables <- c(colnames(soil_variable))

categoris_and_databases <- hits_and_properties %>%
  select(category, database_clean) %>%
  distinct(database_clean, .keep_all = TRUE)

# 1. Compute correlations and p-values
safe_cor <- function(x, y) {
  if (sum(is.finite(x) & is.finite(y)) >= 3) cor(x, y, method = "spearman", use = "pairwise.complete.obs") else NA
}

safe_p <- function(x, y) {
  if (sum(is.finite(x) & is.finite(y)) >= 3) cor.test(x, y, method = "spearman")$p.value else NA
}

cor_results <- hits_and_properties %>%
  group_by(database_clean) %>%
  dplyr::summarise(
    dplyr::across(
      all_of(soil_variables),
      list(
        cor = ~ safe_cor(count, .),
        p = ~ safe_p(count, .)
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  ) %>%
  left_join(
    categoris_and_databases,
    by = "database_clean"
  )


# 2. Prepare matrices

# Separate correlation and p-value matrices
cor_matrix <- cor_results %>%
  mutate(row_id = paste(database_clean, sep = "_")) %>%
  select(row_id, starts_with("cor_")) %>%
  tibble::column_to_rownames("row_id") %>%
  as.matrix()

colnames(cor_matrix) <- gsub("^cor_", "", colnames(cor_matrix))

p_matrix <- cor_results %>%
  mutate(row_id = paste(database_clean, sep = "_")) %>%
  select(row_id, starts_with("p_")) %>%
  tibble::column_to_rownames("row_id") %>%
  as.matrix()

# Replace NAs
cor_matrix[is.na(cor_matrix)] <- 0
p_matrix[is.na(p_matrix)] <- 1

# 3. Collapsed significance per taxon
sig_labels_taxa <- cor_results %>%
  group_by(database_clean) %>%
  dplyr::summarise(across(starts_with("p_"), ~ min(.x, na.rm = TRUE))) %>%
  column_to_rownames("database_clean") %>%
  as.matrix()

sig_labels_taxa <- ifelse(sig_labels_taxa <= 0.001, "***",
                          ifelse(sig_labels_taxa <= 0.01, "**",
                                 ifelse(sig_labels_taxa <= 0.05, "*", "")))

# 4. Category annotation
categories <- unique(cor_results$category)

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
  Category = cor_results$category,
  col = list(Category = category_colors),
  show_annotation_name = FALSE
)


# 7. Plot heatmap with significance
col_fun <- colorRamp2(c(-1, 0, 1), c("#D6604D", "white", "#4993C3"))

norway_gene_corr <- Heatmap(
  cor_matrix,
  name = "Spearman's rho",
  col = col_fun,
  na_col = "grey90",
  show_row_names = T,
  show_column_names = TRUE,
  left_annotation = row_ha,
  cluster_rows = TRUE,
  #row_split = cor_results$category,
  row_title = NULL,
  column_names_side = "bottom",
  row_names_max_width = unit(75, "mm"),
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, col = "black"), # this is to avoid overlapping with labels
  column_names_gp = gpar(fontsize = 8),
  layer_fun = function(j, i, x, y, w, h, fill) { for(idx in seq_along(i)) {
    # Extract taxon ignoring P.level
    taxon <- rownames(cor_matrix)[i[idx]]
    # Pull significance from collapsed matrix
    if(taxon %in% rownames(sig_labels_taxa)) { grid.text(sig_labels_taxa[taxon, j[idx]], x[idx], y[idx], gp = gpar(fontsize = 8)) } } } )
