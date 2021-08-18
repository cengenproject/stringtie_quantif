library(GenomicFeatures)
library(tidyverse)
library(wbData)

options(wb_dir_cache = "/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277/")

gids <- wb_load_gene_ids(277)
txdb <- wb_load_TxDb(277)


# On the result of collapse_tx.R, load the GRanges and make a GFF and TxDb out of it.

new_exons_db <- readRDS("intermediates/210819_new_db_as_granges.rds")


#~ Build metadata columns ----
mcols(new_exons_db)$Name <- mcols(new_exons_db)$exon_name

mcols(new_exons_db)$transcript_name <- mcols(new_exons_db)$exon_name %>% 
  stringr::str_match("^([A-Za-z0-9.+_]+)\\.e[0-9]+$") %>%
  (function(x) x[,2])()


roworder <- mcols(new_exons_db)$transcript_name %>%
  stringr::str_match("^([A-Za-z0-9_]+\\.t{0,1}[0-9]+)[a-z.+0-9]*$") %>%
  magrittr::extract(,2) %>%
  match(gids$sequence)


mcols(new_exons_db)$gene_id <- gids$gene_id[roworder]
mcols(new_exons_db)$gene_symbol <- gids$symbol[roworder]








# Create ids for the new exons and transcripts
exid <- mcols(new_exons_db)$exon_id
exid[is.na(exid)] <- seq(from = 1+max(exid, na.rm = TRUE),
                         to = (max(exid, na.rm = TRUE) + length(which(is.na(exid)))))
mcols(new_exons_db)$exon_id <- exid

mcols(new_exons_db)$transcript_id <- as.integer(as.factor(mcols(new_exons_db)$transcript_name))

# finalize metadata columns
mcols(new_exons_db)$type <- 'exon'

mcols(new_exons_db)$source <- 'WS277_collapsed'
genome(new_exons_db) <- "WS277_collapsed"


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


#~ Export ----

saveDb(col_txdb, "intermediates/210818_collapsed_WS277.txdb.sqlite")
rtracklayer::export(transcriptsBy(col_txdb), "intermediates/210818_collapsed_annotation_WS277.gff3")
