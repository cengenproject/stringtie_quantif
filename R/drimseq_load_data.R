# DRIMSeq preprocessing
# Updated to use stc_q from 220322



## Initializations ----

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(tidyverse)
  library(DRIMSeq)
})


## Prepare Data ----

# From StringTie2 (stc_q)
path_data <- "intermediates/220322_str_q_outs/summaries"

# transcripts table
transcripts_tbl <- readr::read_tsv(file.path(path_data, "transcripts_table.tsv"),
                                   col_types = "cddcccc")



# sample counts
cnts <- readr::read_csv(file.path(path_data, "transcript_count_matrix.csv"),
                        show_col_types = FALSE) |>
  add_column(gene_id = transcripts_tbl$gene_id[match(cnts$transcript_id,
                                                     transcripts_tbl$transcript_id)],
             .after = "transcript_id") |>
  dplyr::rename(feature_id = transcript_id)



# samples table
samples_table <- colnames(cnts) |>
  setdiff(c("feature_id", "gene_id")) |>
  enframe(value = "sample_id", name = NULL) |>
  mutate(neuron_id = str_match(sample_id,
                               pattern = "^([A-Z0-9]{1,4}|Ref)r[0-9]{1,4}$")[,2])

any(is.na(samples_table$neuron_id))




# Create DRIMSeq object
dms <- dmDSdata(as.data.frame(cnts),
                as.data.frame(samples_table))

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
# after filtering,  7,237 / 46,925 genes left
#                  21,581 / 64,758 transcripts left


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

qs::qsave(fitdms, "intermediates/2023-03-30_drimseq_fitdms.rds")

