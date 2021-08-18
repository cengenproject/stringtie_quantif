# Idea: collapse transcripts that "look" similar

# Principle:
#  Compare all tx 2 by 2. If high overlap, combine them


## Initializations ----
library(GenomicFeatures)
library(tidyverse)
library(wbData)

options(wb_dir_cache = "/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277/")

gids <- wb_load_gene_ids(277)
txdb <- wb_load_TxDb(277)



# Functions ----
get_equivalent <- function(tx1, tx2,
                           threshold_tot_prop = 0.05,
                           threshold_tot_bp = 20,
                           threshold_exon_prop = 0.20,
                           threshold_exon_bp = 20){
  # Consider these transcripts equivalent if the total difference is less than treshold_tot
  # and for every exon, the difference is less than treshold_exon.
  # For both, indicating a threshold in proportion as well as bp, one is enough.
  # If equivalent return the combined GRangesList.
  
  if(NROW(tx1) != NROW(tx2)){
    # Finished: don't have the same exons
    new_set <- GRangesList(tx1, tx2)
  } else if(! identical(findOverlaps(tx1, tx2, select="first"), 1:NROW(tx1))){
    # Finished: don't have the same exons
    new_set <- GRangesList(tx1, tx2)
  } else{
    # get proportion of overlap between each exon pair
    # Note: can't use pintersect on GRanges when one exon is
    # contained within another (e.g. WBGene00000041).
    # Use GRangeLists as in https://support.bioconductor.org/p/64161/
    exon_bp <- tryCatch(width(psetdiff(tx1,tx2)),
                        error= function(e) NA)
    if(length(exon_bp) == 1 && is.na(exon_bp)){
      exon_bp <- width(BiocGenerics::setdiff(as(tx1, "GRangesList"),
                               as(tx2, "GRangesList")))
    }
    
    
    exon_prop <- exon_bp/width(punion(tx1, tx2))
    
    
    # Get proportion of overlap on total transcripts lengths
    tot_bp <- sum(width(BiocGenerics::setdiff(tx1, tx2)))
    tot_prop <- tot_bp / sum(width(BiocGenerics::union(tx1, tx2)))
    
    
    if((tot_prop < threshold_tot_prop || tot_bp < threshold_tot_bp) &&
       (max(exon_prop) < threshold_exon_prop || max(exon_bp) < threshold_exon_bp)){
      # Finished: combine equivalent transcripts
      new_set <- GRangesList(BiocGenerics::union(tx1, tx2))
      
      # Make names
      seq_name <- paste0(stringr::str_match(mcols(tx1)$exon_name[1],
                                            "^([A-Za-z0-9_]+\\.[0-9]{1,2})[a-z+.0-9]+e[0-9]+")[2],
                         stringr::str_match(mcols(tx1)$exon_name[1],
                                            "^[A-Za-z0-9_]+\\.[0-9]{1,2}([a-z+.0-9]+)e[0-9]+")[2],
                         "+",
                         stringr::str_match(mcols(tx2)$exon_name[1],
                                            "^[A-Za-z0-9_]+\\.[0-9]{1,2}([a-z+.0-9]+)e[0-9]+")[2]
      )
      # Since we can have isoforms (letters) or only transcripts, remove the useless dots
      seq_name <- stringr::str_replace_all(seq_name,
                                           c("\\.\\." = "\\.", "\\.\\+" = "\\+", "\\+\\." = "\\+", "\\.$"=""))
      
      mcols(new_set[[1]])$exon_name <- paste0(seq_name, ".e", seq(NROW(new_set[[1]])))
      mcols(new_set[[1]])$exon_rank <- seq(NROW(new_set[[1]]))
    } else{
      # Finished: not equivalent bc prop above threshold
      new_set <- GRangesList(tx1, tx2)
    }
  }
  return(new_set)
}



simplify_gene <- function(gene){
  tx_set <- txdf$TXID[txdf$GENEID == gene]
  if(length(tx_set) > 1){
    # initialize with first 2 transcripts
    new_ex_set <- get_equivalent(all_exons[[ tx_set[1] ]], all_exons[[ tx_set[2] ]])
    # Compare all others to existing
    next_tx <- 3
    while(next_tx <= length(tx_set)){
      offset <- 0
      nb_existing <- length(new_ex_set)
      for(existing in 1:nb_existing){
        upd <- get_equivalent(new_ex_set[[existing-offset]], all_exons[[ tx_set[next_tx] ]])
        if(length(upd) == 1){
          # Combined with an existing
          new_ex_set[[existing-offset]] <- NULL
          offset <- offset + 1 # need to offset, as we just deleted one element
          new_ex_set <- c(new_ex_set, upd)
          break
        } else if(existing == nb_existing){
          # No existing equivalent found, add as is
          new_ex_set <- c(new_ex_set, GRangesList(all_exons[[ tx_set[next_tx] ]]))
        }
      }
      next_tx <- next_tx +1
    }
  } else{
    # Single transcript: keep it
    new_ex_set <- GRangesList(all_exons[[ tx_set[1] ]])
  }
  return(new_ex_set)
}





# Format data ----
txdf <- AnnotationDbi::select(txdb,
                              keys=keys(txdb, keytype = "GENEID"),
                              columns=c("GENEID","TXNAME", "TXID"),
                              keytype = "GENEID")

genelist <- keys(txdb, "GENEID")
all_exons <- exonsBy(txdb, by="tx")


# Do for each gene ----

cat("Starting run, may take 40 min\n\n")


# # To run on laptop
# new_exons_db <- lapply(genelist[1:3], simplify_gene)

# To run on SLURM
new_exons_db <- parallel::mclapply(genelist, simplify_gene, mc.cores = 15)


new_exons_db <- unlist(unlist(List(new_exons_db)))

saveRDS(new_exons_db, "intermediates/210819_new_db_as_granges.rds")
cat("Saved\n")


#~ Build metadata columns ----
mcols(new_exons_db)$transcript_name <- mcols(new_exons_db)$exon_name |> 
  stringr::str_match("^([A-Za-z0-9.+_]+)\\.e[0-9]+$") |>
  (\(x) x[,2])()


mcols(new_exons_db)$gene_id <- mcols(new_exons_db)$transcript_name |>
  stringr::str_match("^([A-Za-z0-9_]+\\.[t0-9]+)[a-z.+0-9]*$") |>
  (\(x) x[,2])() |>
  wb_seq2id(gids, warn_missing = TRUE)

mcols(new_exons_db)$gene_symbol <- mcols(new_exons_db)$transcript_name |> 
  stringr::str_match("^([A-Z0-9]+\\.[t0-9]+)[a-z.+0-9]*$") |> 
  (\(x) x[,2])() |>
  wb_seq2name(gids, warn = TRUE)


col_txdb <- makeTxDbFromGRanges(new_exons_db)



#~ Checks ----
# nb transcripts
cat("Number of transcripts: ")
mcols(new_exons_db)$transcript_name |> 
  unique() |> 
  length() |> 
  cat()
cat("\nNumber of genes: ")
mcols(new_exons_db)$gene_id |> 
  unique() |> 
  length()
cat("\n vs initial number of genes: ",length(genelist),"\n")


#~ Export ----

saveDb(col_txdb, "intermediates/210818_collapsed_WS277.txdb.sqlite")
rtracklayer::export(transcriptsBy(col_txdb), "intermediates/210818_collapsed_annotation_WS277.gff3")

