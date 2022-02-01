# Read the StringTie output, make summary files with TPM etc
# to be called on the cluster

# call with arguments: first, quantification dir, second output dir

# Init ----

library(tidyverse)


# qdir <- "intermediates/220128_str_q_outs/quantifications"
# expdir <- "intermediates/220128_str_q_outs/summaries"

args <- commandArgs(TRUE)
if(length(args) != 2L){
  stop("Incorrect number of arguments, expected 2, got ", length(args))
}

qdir <- args[[1]]
expdir <- args[[2]]

# check arguments
if(!dir.exists(qdir)){
  stop("Non-existent quantifications directory: ", qdir)
}
if(!dir.exists(expdir)){
  warning("Creating non-existent output directory: ", expdir)
  dir.create(expdir)
}

if(length(list.dirs(qdir)) < 1L){
  stop("Quantifications dir empty: ", qdir)
}

if(any(file.exists(file.path(expdir, "samples_table.tsv"),
                   file.path(expdir, "transcripts_table.tsv"),
                   file.path(expdir, "tx_cov.tsv"),
                   file.path(expdir, "tx_FPKM.tsv"),
                   file.path(expdir, "tx_TPM.tsv")))){
  warning("Will overwrite existing output files")
  Sys.sleep(3)
}




# Read in data ----
samples_table <- tibble(sample_id = list.dirs(qdir,full.names=FALSE,recursive = FALSE),
                        sample_path = list.dirs(qdir,full.names=TRUE,recursive = FALSE))

stopifnot(all(! is.na(str_match(samples_table$sample_id, "^([A-Z0-9ef]{1,4}r[0-9]{1,3})$")[,2])))



import_file <- function(sample_id,sample_path){
  imp <- rtracklayer::import(paste0(sample_path, "/", sample_id, ".gtf"),
                             feature.type = "transcript") %>%
    as_tibble() %>%
    mutate(across(c(cov,FPKM,TPM), as.double))
  
  if("gene_name" %in% colnames(imp)){
    # I5r210 has no gene_name column, somehow
    imp <- imp %>%
      mutate(gene_name = if_else(is.na(ref_gene_name), gene_name, ref_gene_name)) %>%
      select(-c(width, source, type, score, phase, ref_gene_name))
  } else{
    imp <- imp %>%
      mutate(gene_name = ref_gene_name) %>%
      select(-c(width, source, type, score, phase, ref_gene_name))
  }
  imp
}


samples_table$content <- purrr::pmap(samples_table,
                                     import_file)

samples_table$sample_path <- NULL


# all samples have the same transcripts (as we used ref annotation)
stopifnot(length(unique(map_int(samples_table$content, nrow))) == 1L)
stopifnot(all(map_int(samples_table$content, ncol) == 10L))




# Create outputs ----
#    * samples table (already done)
#    * transcripts table
#    * cov matrix
#    * FPKM matrix
#    * TPM matrix



transcripts_table <- samples_table$content[[1]] %>%
  select(-c(cov,FPKM,TPM))



vals_long <- pmap_dfr(samples_table,
         function(sample_id, content){
           content %>%
             select(transcript_id, cov, FPKM, TPM) %>%
             add_column(sample_id = sample_id)
           })


mat_cov <- vals_long %>%
  select(-FPKM,-TPM) %>%
  pivot_wider(transcript_id,
              names_from = sample_id,
              values_from = cov)

mat_fpkm <- vals_long %>%
  select(-cov,-TPM) %>%
  pivot_wider(transcript_id,
              names_from = sample_id,
              values_from = FPKM)

mat_tpm <- vals_long %>%
  select(-FPKM,-cov) %>%
  pivot_wider(transcript_id,
              names_from = sample_id,
              values_from = TPM)


samples_table <- samples_table %>%
  select(-content) %>%
  mutate(neuron_id = str_match(sample_id, "^([A-Z0-9ef]{1,4})r[0-9]{1,3}$")[,2])



# Write files ----
write_tsv(transcripts_table,
          file.path(expdir, "transcripts_table.tsv"))

write_tsv(samples_table,
          file.path(expdir, "samples_table.tsv"))


write_tsv(vals_long,
          file.path(expdir, "tx_long.tsv"))


write_tsv(mat_cov,
          file.path(expdir, "tx_cov.tsv"))

write_tsv(mat_fpkm,
          file.path(expdir, "tx_FPKM.tsv"))

write_tsv(mat_tpm,
          file.path(expdir, "tx_TPM.tsv"))



