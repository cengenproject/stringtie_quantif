suppressPackageStartupMessages(library(tidyverse))

library(wbData)
gids <- wb_load_gene_ids(277)

source("R/utils.R")



# Do 2 things:
# * from the collapsed GTF, find gene name for each annotated tx
# * from the merged GTF, attribute novels tx into genes (novel or existing)


# 1/  gene names from collapsed ----
#Use the GFF created by collapsing to make a tx <-> gene look-up table.


## Load separately tx and exons (bc tx entries have gene ID, exons entries have descriptive exon name
# created during the collapsing)
coll_tx <- rtracklayer::readGFF("data/intermediates_collapse/210909_collapsed_annotation_WS277.gff3",
                                filter = list(type="mRNA"),
                                tags = c("ID","Parent")) |>
  as_tibble() |>
  mutate(gene_id = str_match(Parent, "^GeneID:(WBGene[0-9]{8})$")[,2]) |>
  select(-Parent, transcript_nb = ID)

# Check that unique Parent for each tx
# map_int(coll_tx$Parent, length) |> table()

coll_exon <- rtracklayer::readGFF("data/intermediates_collapse/210909_collapsed_annotation_WS277.gff3",
                                  filter = list(type="exon"),
                                  tags = c("Name","Parent")) |>
  as_tibble() |>
  mutate(tx_name = stringr::str_match(Name,"^([A-Za-z0-9.+_]+)\\.e[0-9]+$")[,2],
         transcript_nb = unlist(Parent)) |>
  select(-Parent)

# Check that unique Parent for each exon
# map_int(coll_exon$Parent, length) |> table()



## Combine them into a single look-up table

tx2g_annot <- full_join(
  select(coll_tx, transcript_nb, gene_id),
  select(coll_exon, transcript_nb, tx_name) |> distinct(),
  by = "transcript_nb")



# 2/ gene grouping of novel tx ----

# From the merged.gtf created by StringTie.
# Contains 2 sources: "rtracklayer" for tx that were not really detected in data, it's kept in the GTF unchanged.
# If tx was detected, gene got a new "MSTRG.x" id (but old gene kept as "ref_gene_id").
# If previously annotated, the tx ID is kept (from collapsed, ie with format "TxID:xxxx"), if novel
# gets a name based on gene ID, eg MSTRG1.1. Note: several novel transcripts can be attached to annotated gene,
# only the annotated transcript keeps its name and ref_gene_id, for the novel ones one has to look at the MSTRG gene_id.


# Load transcripts from merge (annotated and novel)
merg_tx <- rtracklayer::readGFF("data/210929_stc_nq_summaries/merged.gtf",
                                filter = list(type="transcript"),
                                tags = c("gene_id","transcript_id","ref_gene_id")) |>
  as_tibble()

# get gene for each tx, clean gene name for annotated ones
# group together annotated genes that have the same ST gene
gene_annot <- merg_tx |>
  select(gene_id, ref_gene_id) |>
  arrange(gene_id, ref_gene_id) |>
  group_by(gene_id) |>
  fill(ref_gene_id) |>
  distinct() |>
  ungroup() |>
  mutate(ref_gene_id = str_match(ref_gene_id, "^GeneID:(WBGene[0-9]{8})$")[,2],
         ref_gene_seq = wb_id2seq(ref_gene_id, gids),
         ref_gene_name = i2s(ref_gene_id, gids)) |>
  group_by(gene_id) |>
  summarize(ref_gene_id = if_else(all(is.na(ref_gene_id)),
                                  NA_character_,
                                  paste0(ref_gene_id, collapse = "+")),
            ref_gene_seq = if_else(all(is.na(ref_gene_id)),
                                    NA_character_,
                                    paste0(ref_gene_seq, collapse = "+")),
            ref_gene_name = if_else(all(is.na(ref_gene_name)),
                                    NA_character_,
                                    paste0(ref_gene_name, collapse = "+")),
            .groups = "drop")

# associate each tx with gene
all_tx <- merg_tx |>
  select(transcript_id, gene_id) |>
  left_join(gene_annot, by = "gene_id")

# give ID to novel genes
id_for_novel_genes <- all_tx |>
  filter(is.na(ref_gene_id)) |> 
  select(gene_id) |>
  distinct() |>
  arrange() |>
  mutate(ref_gene_id = create_gene_id(str_match(gene_id, "^MSTRG\\.([0-9]{1,5})$")[,2]),
         ref_gene_seq = create_seq_name(str_match(gene_id, "^MSTRG\\.([0-9]{1,5})$")[,2])) |>
  column_to_rownames("gene_id")



# Update table with all gene names (both annotated and novel)
all_tx <- all_tx |>
  mutate(ref_gene_id = if_else(is.na(ref_gene_id),
                                id_for_novel_genes[gene_id,"ref_gene_id"],
                                ref_gene_id),
         ref_gene_seq = if_else(is.na(ref_gene_seq),
                               id_for_novel_genes[gene_id,"ref_gene_seq"],
                               ref_gene_seq))



## make tx table

# init with tx and gene names
tx_table <- merg_tx |>
  select(transcript_id, gene_id) |>
  left_join(all_tx, by = "transcript_id")

# check
table(tx_table$gene_id.x == tx_table$gene_id.y)
tx_table <- tx_table |>
  select(-gene_id.x,
         gene_id = gene_id.y)

# add tx name when known
tx_table <- left_join(tx_table,
                      tx2g_annot |>
                        select(-gene_id),
                      by = c(transcript_id = "transcript_nb"))

# check that tx_name is gene seq name + something at the end
xx <- tx_table |> filter(!is.na(tx_name),
                         !str_detect(ref_gene_id, "\\+"))
table(xx$ref_gene_seq != str_match(xx$tx_name, "^([A-Z0-9cel_]+\\.t?[0-9]{1,4})")[,2])
rm(xx)


# create tx name for novel ones
tx_table <- tx_table |>
  group_by(gene_id) |>
  mutate(tx_name = if_else(is.na(tx_name),
                           paste0(ref_gene_seq, ".s", row_number()),
                           tx_name))


# 3/ export ----


write_csv(tx_table,
          "data/210929_stc_nq_summaries/2021-10-14_tx_table.csv")

