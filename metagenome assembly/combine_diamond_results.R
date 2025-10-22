library(dplyr)
library(readr)
library(purrr)

# Path to all samples
samples <- list.dirs("/home/bogki545/00_nesi_projects/uoo04534_nobackup/Norway/3_assembled", recursive = FALSE, full.names = TRUE)


# Loop over samples
all_hits <- map_dfr(samples, function(samp_dir) {
  sample_name <- basename(samp_dir)
  hits_files <- list.files(file.path(samp_dir, "greening_hits"), pattern="\\.tsv$", full.names=TRUE)
  
  map_dfr(hits_files, function(f) {
    db_name <- gsub("_hits\\.tsv$", "", basename(f))
    read_tsv(f, col_names = c("gene_id", "hit_id", "pident", "length", "evalue", "bitscore"),
             col_types = cols(
               gene_id = col_character(),
               hit_id = col_character(),
               pident = col_double(),
               length = col_double(),
               evalue = col_double(),
               bitscore = col_double()
             )) %>%
      mutate(
        sample = sample_name,
        database = db_name
      )
  })
})

saveRDS(all_hits, file = "/home/bogki545/00_nesi_projects/uoo04534_nobackup/Norway/greening_results/all_hits_new.rds")
