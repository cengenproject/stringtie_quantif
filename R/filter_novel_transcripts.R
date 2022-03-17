# Clean up stringtie discovery

# First, run str_sc_n using mix of short and long reads to discover new transcripts.
# But many of those tx are "fusion" from different genes, probably not true
# Here, load merged and:
#   * Identify reference tx
#   * Identify novel tx which have only small change compared to ref tx
#   * Identify wrong tx which are fusion of two ref tx with a bit more
#   * Identify fully novel tx from novel genes (e.g. STRG.1)
#   * Add back ref tx not found in merged (incl pseudogenes, ncRNA, ...)

suppressPackageStartupMessages(library(GenomicFeatures))
library(wbData)

gr_mgd <- rtracklayer::import("data/2022-03-14_str_sc_n_merged.gtf")

gr_novel <- gr_mgd[strand(gr_mgd) != "*" &
                     is.na(gr_mgd$reference_id), ]


# prep data
novel_exons <- exonsBy(makeTxDbFromGRanges(gr_novel), use.names = TRUE)

ref_exons <- exonsBy(wbData::wb_load_TxDb(281), use.names = TRUE)


# prop of the exonic seq that overlaps with ref tx
mgd <- mergeByOverlaps(novel_exons, ref_exons)
prop_overlap <- sum(width(intersect(mgd$novel_exons, mgd$ref_exons)))/sum(width(mgd$novel_exons))

nb_ov <- countOverlaps(novel_exons, ref_exons)



## Novel genes ----
# if less than 10% of its exonic seq overlaps with sum ref tx

prop_single_overlap <- prop_overlap[names(nb_ov)[nb_ov == 1]]


tx_from_novel_gene <- names(nb_ov)[nb_ov == 0] |>
  union(names(prop_single_overlap[prop_single_overlap < .1]))


## Alternative tx ----
# if overlaps with set of ref tx from single gene, by more than 30% of its exonic seq


tx2g_tab <- wb_load_tx2gene(281)

new_tx_match_g <- data.frame(new_tx = names(mgd$novel_exons),
                             genes = wb_tx2g(names(mgd$ref_exons), tx2g_tab, warn_missing = TRUE))

nb_separate_genes_overlapping <- aggregate(genes ~ new_tx,
                                           data = new_tx_match_g[prop_overlap > 0.1,],
                                           FUN = \(.x) length(unique(.x)))

tx_alt_existing_gene <- nb_separate_genes_overlapping$new_tx[nb_separate_genes_overlapping$genes == 1] |>
  intersect(names(prop_overlap)[prop_overlap > .3])


## Fusion btw genes ----
# if overlaps with 2 or more ref tx from different genes, and contains more than 30% of the exonic seq of at least 1 tx


ref_prop_overlap <- sum(width(intersect(mgd$novel_exons, mgd$ref_exons)))/sum(width(mgd$ref_exons))


nb_separate_genes_overlapping_ref50 <- aggregate(genes ~ new_tx,
                                                 data = new_tx_match_g[ref_prop_overlap > 0.3,],
                                                 FUN = \(.x) length(unique(.x)))


tx_fusion <- nb_separate_genes_overlapping_ref50$new_tx[nb_separate_genes_overlapping_ref50$genes > 1]



## Finally ----

tx_alt_existing_gene <- tx_alt_existing_gene |>
  setdiff(tx_fusion)

# Find parent gene

tx2g_new_tx <- new_tx_match_g[prop_overlap > 0.1 &
                                new_tx_match_g$new_tx %in% tx_alt_existing_gene,] |>
  dplyr::group_by(new_tx) |>
  dplyr::summarize(gene = unique(genes))



### Classify examples ----
yy <- character(length(xx))

yy[xx %in% tx_alt_existing_gene] <- "alt"
yy[xx %in% tx_fusion] <- "fusion"

yy[xx %in% tx_from_novel_gene] <- "novel"
clipr::write_clip(yy)


## Export ----

# make GTF manually (as rtracklayer collapses exons from different tx)
txdf <- gr_novel[gr_novel$type == "transcript" &
                   gr_novel$transcript_id %in% tx_alt_existing_gene] |>
  as.data.frame() |>
  dplyr::mutate(source = "StringTie",
                type = "transcript",
                score = "1000",
                frame = ".",
                gene_id = tx2g_new_tx$gene[match(transcript_id, tx2g_new_tx$new_tx)],
                attributes = paste0('gene_id "',gene_id,'"; transcript_id "',
                                    transcript_id,'";')) |>
  dplyr::select(chrom = seqnames,
                source,
                type,
                start,
                end,
                score,
                strand,
                frame,
                attributes)

exonsdf <- gr_novel[gr_novel$type == "exon" &
                          gr_novel$transcript_id %in% tx_alt_existing_gene] |>
  as.data.frame() |>
  dplyr::group_by(transcript_id) |>
  dplyr::mutate(source = "StringTie",
                type = "exon",
                score = "1000",
                frame = ".",
                gene_id = tx2g_new_tx$gene[match(transcript_id, tx2g_new_tx$new_tx)],
                attributes = paste0('gene_id \"',gene_id,'\"; transcript_id \"',
                                    transcript_id,'\"; ',
                                    'exon_number \"',dplyr::row_number(),'\";')) |>
  dplyr::ungroup() |>
  dplyr::select(chrom = seqnames,
                source,
                type,
                start,
                end,
                score,
                strand,
                frame,
                attributes)


dplyr::bind_rows(txdf, exonsdf) |>
  readr::write_tsv("data/220317_novel.gtf",
                   col_names = FALSE,
                   escape = "none")


# check result

shell("cat data/220317_novel.gtf | cut -f3 | sort | uniq -c")

length(unlist(novel_exons[tx_alt_existing_gene]))
length(tx_alt_existing_gene)


# merge back the reference transcripts
shell(paste0("cat ", wb_get_gtf_path(281) |>
               stringr::str_remove(".gz")," >> data/220317_novel.gtf"))

shell("sort -k1,1 -k4,4n data/220317_novel.gtf > data/220317_novel.sorted.gtf")


