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
path_raw <- "raw/stc_dcq3/"
path_data <- "data/stc_dcq3_v4"




# get tx <-> gene table from "find_gene_for_tx4.R": need to run separately
transcripts_tbl <- readRDS("data/transcript_WBGene_lut.rds")


# Add the new unnamed genes to the neuron expression matrix, obtained by running "gene_calling.R"
ungenes <- transcripts_tbl$gene_id[!is.na(str_extract(transcripts_tbl$gene_id, "UNG"))]

expr_gene_by_neuron <- readRDS("data/expression_22samples_thres14.rds")

expr_gene_by_neuron <- rbind(expr_gene_by_neuron,
                             matrix(nrow = length(ungenes),
                                    ncol = ncol(expr_gene_by_neuron),
                                    dimnames = list(ungenes,
                                                    colnames(expr_gene_by_neuron))))

# Data from STAR alignment + collapsedGTF + StringTie2 (bss1cst1_1 pipeline)
# Read sample counts
cnts <- readr::read_csv(file.path(path_raw, "transcript_count_matrix.csv")) %>%
  as.data.frame()

# rename genes and transcripts
cnts$gene_id <- transcripts_tbl$gene_id[match(cnts$transcript_id, transcripts_tbl$feature_id)]
cnts <- rename(cnts, feature_id = transcript_id)


# Prepare samples table
old_samples_table <- readr::read_tsv("raw/full_sample_list.tsv", col_names = c("sample_id", "neuron"))%>%
  mutate(col_name = paste0("X",gsub("-", ".", sample_id))) %>%
  as.data.frame()

new_sample_table <- str_match(colnames(cnts), "^([A-Z0-9]{1,4}|Ref)r[0-9]{1,4}$") %>%
  as_tibble() %>%
  select(sample_id = V1,
         neuron = V2,
         col_name = V1) %>%
  filter(!is.na(sample_id))

samples_table <- rbind(old_samples_table, new_sample_table)


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
# after filtering,  3,700 / 20,000 genes left
#                  9,000 / 45,000 transcripts left


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

saveRDS(fitdms, "data/2020-04-20_fitdms")

