library(dplyr)
library(tidyr)
library(stringr)

# import DIAMNOD results
greening_hits <- readRDS("~/.rds")
# set plot IDs
plotIDs <- c("6109_072_S129" = "A-K",
             "6109_073_S130" = "G-K",
             "6109_074_S131" = "H-K",
             "6109_075_S132" = "E-K",
             "6109_076_S133" = "B-K")
# rename genes
gene_renames <- c(
  "Acetyl CoA synthase AcsB \\(99\\)" = "Acetyl CoA synthase (AcsB)",
  "Rubisco RbcL 95 \\(Filtered\\)" = "Rubisco (RbcL)",
  "Nitrogenase NifH \\(99\\)" = "Nitrogenase (NifH)",
  "Cytochrome cd1 nitrite reductase NirS \\(sequence\\)" = "Cytochrome cd1 nitrite reductase (NirS)",
  "Isoprene monooxygenase subunit \\(A\\)" = "Isoprene monooxygenase alpha-subunit (isoA)"
)

# filter, assign plot IDs, rename genes, cluster them
greening_hits_filt <- greening_hits %>%
  filter(pident >= 60, length >= 40, evalue <= 1e-5) %>%
  mutate(sample = str_replace_all(sample, plotIDs),
         sample = factor(sample, levels = c("A-K", "B-K", "E-K", "G-K", "H-K"))) %>%
  mutate(database_clean = database %>%
           str_replace_all("_", " ") %>%
           str_replace(regex(" sequences$", ignore_case = TRUE), "") %>%
           str_replace("(.*) ([^ ]+)$", "\\1 (\\2)") %>%
           str_replace_all(gene_renames) %>%
           str_trim()) %>%
  mutate(category = case_when(
    # aerobic respiration
    database_clean %in% c("Cytochrome aa3 oxidase (CoxA)",
                          "Cytochrome bd oxidase (CydA)",
                          "Succinate dehydrogenase (SdhA)",
                          "F-type ATP synthase (AtpA)",
                          "Fumarate reductase (FrdA)",
                          "NADH-ubiquinone oxidoreductase (NuoF)",
                          "Cytochrome cbb3 oxidase (CcoN)",
                          "Cytochrome bo3 oxidase (CyoA)") ~ "Aerobic respiration",
    # alternative e acceptor
    database_clean %in% c("Arsenite oxidase (ARO)",
                          "Decaheme iron reductase (MtrB)",
                          "Reductive dehalogenase (RdhA)",
                          "Iron oxidising cytochrome (Cyc2)") ~ "Alternative e acceptor",
    # alternative e donor
    database_clean %in% c("Arsenate reductase (ArsC)",
                          "Formate dehydrogenase (FdhA)",
                          "Selenate reductase (YgfK)") ~ "Alternative e donor",
    # carbon fixation
    database_clean %in% c("ATP-citrate lyase (AclB)",
                          "Rubisco (RbcL)",
                          "Thaumarchaeota 4-hydroxybutyryl-CoA synthetase (HbsT)",
                          "Acetyl CoA synthase (AcsB)",
                          "Malonyl-CoA reductase (Mcr)") ~ "Carbon fixation",
    # nitrogen metabolism
    database_clean %in% c("Ammonia monooxygenase (amoA)",
                          "Ammonia-forming nitrite reductase (NrfA)",
                          "Copper containing nitrite reductase (NirK)",
                          "Dissimilatory nitrate reductase (NarG)",
                          "Periplasmic nitrate reductase (NapA)",
                          "Nitrite oxidoreductase (NxrA)",
                          "Nitrous oxide reductase (NosZ)",
                          "Ammonia-forming nitrite reductase (NrfA)",
                          "Cytochrome cd1 nitrite reductase (NirS)",
                          "Nitric oxide reductase (NorB)",
                          "Nitrogenase (NifH)") ~ "Nitrogen metabolism",
    # trace gas metabolism
    database_clean %in% c("Carbon monoxide dehydrogenase (CoxL)",
                          "FeFe (hydrogenase)",
                          "NiFe (hydrogenase)") ~ "Trace gas metabolism",
    # phototrophy
    database_clean %in% c("Microbial rhodopsin (RHO)",
                          "Photosystem I reaction centre (PsaA)",
                          "Photosystem II reaction centre (PsbA)") ~ "Phototrophy",
    # sulfur metabolism
    database_clean %in% c("Flavocytochrome c sulfide dehydrogenase (FCC)",
                          "Sulfide quinone oxidoreductase (Sqr)",
                          "Thiosulfohydrolase (SoxB)") ~ "Sulfur metabolism",
    # other
    database_clean %in% c("Isoprene monooxygenase alpha-subunit (isoA)") ~ "Other"
  )) %>%
  mutate(category = factor(category,
                           levels = c("Aerobic respiration",
                                      "Alternative e acceptor",
                                      "Alternative e donor",
                                      "Carbon fixation",
                                      "Nitrogen metabolism",
                                      "Sulfur metabolism",
                                      "Trace gas metabolism",
                                      "Phototrophy",
                                      "Other"))) %>%
  mutate(original_n_reads = case_when( # original number of reads
    sample == "A-K" ~ 75760010, # old 55178250
    sample == "G-K" ~ 76954132,
    sample == "H-K" ~ 71900346,
    sample == "E-K" ~ 71605938,
    sample == "B-K" ~ 83080025 # old 56645966
  ))

