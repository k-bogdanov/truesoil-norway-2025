# theme
theme_norway <- theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(color = "black", family = "roboto"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.text = element_text(colour = "black", size = 11, margin = margin(t = 5, r = 5, b = 5, l = 5)),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

set.seed(12358)

# import ps objects
Norway_16S_ps <- readRDS("C:/Users/bogda/Documents/PhD/TRUESOIL/molecular_data/2025/phyloseqs/ps_Norway_16S.rds") %>%
  name_na_taxa(na_label = "Unidentified <tax> (<rank>)")

Norway_ITS_ps <- readRDS("C:/Users/bogda/Documents/PhD/TRUESOIL/molecular_data/2025/phyloseqs/ps_Norway_ITS.rds") %>%
  name_na_taxa(na_label = "Unidentified <tax> (<rank>)")


# remove "p__" etc from fungal names
tax_table(Norway_ITS_ps)[, colnames(tax_table(Norway_ITS_ps))] <- gsub(tax_table(Norway_ITS_ps)[, colnames(tax_table(Norway_ITS_ps))], pattern = "[a-z]__", replacement = "")

# rarefy
Norway_16S_ps_r <- rarefy_even_depth(Norway_16S_ps, 5000, rngseed = TRUE)
Norway_ITS_ps_r <- rarefy_even_depth(Norway_ITS_ps, 5000, rngseed = TRUE)

# pull metadata
Norway_metadata <- data.frame(sample_data(Norway_16S_ps_r))
