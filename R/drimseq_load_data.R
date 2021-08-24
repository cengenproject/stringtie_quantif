# DRIMSeq preprocessing
# Version "stc_dcq3_v4": using "stc_dcq3" data, 4th version of the scripts

# First need to run STAR alignment + collapsedGTF + StringTie2 (bss1cst1_1 pipeline)
# then need to run find_gene_for_tx4.R


## Initializations ----

library(AnnotationDbi)
library(GenomicFeatures)
library(tidyverse)
library(DRIMSeq)


## Prepare Data ----
# From StringTie2 (stc_dcq5)
path_data <- "data/str2_summaries"


# get tx <-> gene table from "find_gene_for_tx5.R": need to run separately
transcripts_tbl <- readRDS("data/transcript_WBGene_lut.rds")




# Data from STAR alignment + collapsedGTF + StringTie2 (bss1cst1_1 pipeline)
# Read sample counts
cnts <- readr::read_csv(file.path(path_data, "transcript_count_matrix.csv")) %>%
  as.data.frame()

# Prepare samples table
old_samples_table <- readr::read_csv("data/210820_full_sample_list.csv", col_types = "ccc") %>%
  dplyr::select(-replicate) %>%
  dplyr::rename(sample_id = sample)

new_sample_table <- str_match(colnames(cnts), "^([A-Z0-9]{1,4}|Ref)r[0-9]{1,4}$") %>%
  as_tibble() %>%
  dplyr::slice(-1) %>%
  dplyr::rename(sample_id = V1,
                neuron = V2)

waldo::compare(old_samples_table%>% arrange(sample_id), new_sample_table%>% arrange(sample_id))

samples_table <- as.data.frame(new_sample_table)




# rename genes and transcripts
cnts$gene_id <- transcripts_tbl$gene_id[match(cnts$transcript_id, transcripts_tbl$feature_id)]
cnts <- dplyr::rename(cnts, feature_id = transcript_id)





# Create DRIMSeq object
dms <- dmDSdata(cnts, samples_table)

# Filter out lowly expressed:
# 1/ we keep only transcripts with at least `min_feature_expr`
#                 in at least `min_samps_feature_expr` samples
# 2/ we keep only transcripts that make up at least `min_feature_prop` of the gene
#                 in at least `min_samps_feature_prop` samples
# 3/ we keep only genes that are expressed at least at `min_gene_expr`
#                 in at least `min_samps_gene_expr` samples
fdms <- dmFilter(dms,
                 min_samps_feature_expr = 3, min_feature_expr = 1,
                 min_samps_feature_prop = 3, min_feature_prop = 0.1,
                 min_samps_gene_expr = 7, min_gene_expr = 20)
# after filtering,  8,791 /  28,697 genes left
#                  49,723 / 113,661 transcripts left


# Prepare contrasts
design_no_int <- model.matrix(~ 0 + neuron, data = DRIMSeq::samples(fdms))

factors_in_design <- colnames(design_no_int)
contr_neurs <- character(0)
for(i in 1:length(factors_in_design)){
  contr_neurs <- c(contr_neurs, paste0(factors_in_design[i],
                                       "-(",
                                       paste(factors_in_design[-i], collapse = "+"),
                                       ")/",
                                       length(factors_in_design)-1))
}


# Fit model ----
pdms <- dmPrecision(fdms, design = design_no_int)
fitdms <- dmFit(pdms, design = design_no_int, verbose = 1)

saveRDS(fitdms, "data/2021-08-20_fitdms.rds")

