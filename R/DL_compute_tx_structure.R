# Export tx structure for DL
# 
# 
# Load GTF. For each tx, keep exons.
# Determine sequence of introns.
# Save one file per gene, with starts and ends of segments and whether exon or intron



## Inits ----


suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(tidyverse)
})

library(wbData)
tx2g_tab <- wb_load_tx2gene(281)


gene_expression <- read.csv("data/thresholded_gene_expression/bsn9_subtracted_integrated_binarized_expression_withVDDD_FDR_0.1_092022.tsv",
                            sep = "\t")
gene_expression2 <- apply(gene_expression, 1, sum)
genes_expressed_in_neurons <- names(gene_expression2)[gene_expression2 > 0]

rm(gene_expression)
rm(gene_expression2)


#~ Make tables of transcript structure ----

txdb <- wb_load_TxDb(281)

exons <- exonsBy(txdb, by = "tx", use.names=TRUE)
introns <- intronsByTranscript(txdb, use.names=TRUE)

all.equal(names(introns), names(exons))


get_tx_structure <- function(tx){
  rbind(exons[[tx]] |>
          as.data.frame() |>
          dplyr::select(chromosome = seqnames,
                        strand = strand,
                        segment_start = start,
                        segment_end = end) |>
          dplyr::mutate(type = "exon"),
        introns[[tx]] |>
          as.data.frame() |>
          dplyr::select(chromosome = seqnames,
                        strand = strand,
                        segment_start = start,
                        segment_end = end) |>
          dplyr::mutate(type = "intron")) |>
    mutate(transcript_id = tx,
           gene_id = wb_tx2g(transcript_id,
                             tx2g_tab,
                             warn_missing = TRUE)) |>
    arrange(segment_start) |>
    relocate(gene_id, transcript_id, .before = chromosome)
}




gene_coords <- wb_load_gene_coords(281) |>
  rename(gene_start = start,
         gene_end = end) |>
  select(-strand, -chr)


segments_coords <- map_dfr(names(exons),
                          get_tx_structure) |>
  left_join(gene_coords,
            by = "gene_id") |>
  mutate(segment_start2 = if_else(strand == "+",
                              as.integer(segment_start - gene_start + 1),
                              as.integer(gene_end - segment_end + 1)),
         segment_end2   = if_else(strand == "+",
                              as.integer(segment_end - gene_start + 1),
                              as.integer(gene_end - segment_start + 1)),
         gene_length = as.integer(gene_end - gene_start +1)) |>
  select(gene_id, transcript_id,
         segment_start = segment_start2,
         segment_end = segment_end2,
         type, gene_length) |>
  arrange(gene_id, transcript_id, segment_start2)

# qs::qsave(segments_coords, "data/intermediates_for_DL/220919_segments_coordinates.qs")
# segments_coords <- qs::qread("data/intermediates_for_DL/220919_segments_coordinates.qs")

export_dir <- "data/intermediates_for_DL/220919_segments_coordinates/"


segments_coords |>
  filter(gene_id %in% genes_expressed_in_neurons) |>
  select(gene_id, transcript_id,
         segment_start = segment_start2,
         segment_end = segment_end2,
         type, gene_length) |>
  split(f = ~gene_id) |>
  iwalk(~write_tsv(.x,
                  paste0(export_dir, .y, ".tsv")))





