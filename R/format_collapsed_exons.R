library(GenomicFeatures)
library(tidyverse)
library(wbData)

options(wb_dir_cache = "/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277/")

gids <- wb_load_gene_ids(277)
txdb <- wb_load_TxDb(277)


# On the result of collapse_tx.R, load the GRanges and make a GFF and TxDb out of it.

new_exons_db <- readRDS("intermediates/210819_new_db_as_granges.rds")


#~ Build metadata columns ----
mcols(new_exons_db)$transcript_name <- mcols(new_exons_db)$exon_name %>% 
  stringr::str_match("^([A-Za-z0-9.+_]+)\\.e[0-9]+$") %>%
  (function(x) x[,2])()


mcols(new_exons_db)$gene_id <- mcols(new_exons_db)$transcript_name %>%
  stringr::str_match("^([A-Za-z0-9_]+\\.[t0-9]+)[a-z.+0-9]*$") %>%
  (function(x) x[,2])() %>%
  wb_seq2id(gids, warn_missing = TRUE)

mcols(new_exons_db)$gene_symbol <- mcols(new_exons_db)$transcript_name %>%
  stringr::str_match("^([A-Z0-9]+\\.[t0-9]+)[a-z.+0-9]*$") %>%
  (function(x) x[,2])() %>%
  wb_seq2name(gids, warn = TRUE)


col_txdb <- makeTxDbFromGRanges(new_exons_db)



#~ Checks ----
# nb transcripts
cat("Number of transcripts: ")
mcols(new_exons_db)$transcript_name %>%
  unique() %>%
  length() %>%
  cat()
cat("\nNumber of genes: ")
mcols(new_exons_db)$gene_id %>%
  unique() %>%
  length()
cat("\n vs initial number of genes: ",length(genelist),"\n")


#~ Export ----

saveDb(col_txdb, "intermediates/210818_collapsed_WS277.txdb.sqlite")
rtracklayer::export(transcriptsBy(col_txdb), "intermediates/210818_collapsed_annotation_WS277.gff3")
