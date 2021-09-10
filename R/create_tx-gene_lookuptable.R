library(tidyverse)

#     Use the GFF created by collapsing to make a tx <-> gene look-up table.



## Load separately tx and exons (bc tx entries have gene ID, exons entries have descriptive exon name
# created during the collapsing)
coll_tx <- rtracklayer::readGFF("data/intermediates_collapse/210909_collapsed_annotation_WS277.gff3",
                                 filter = list(type="mRNA"),
                                 tags = c("ID","Parent")) |>
  as_tibble()

# Check that unique Parent for each tx
map_int(coll_tx$Parent, length) |> table()

coll_exon <- rtracklayer::readGFF("data/intermediates_collapse/210909_collapsed_annotation_WS277.gff3",
                                  filter = list(type="exon"),
                                  tags = c("Name","Parent")) |>
  as_tibble()

# Check that unique Parent for each exon
map_int(coll_exon$Parent, length) |> table()



## Combine them into a single look-up table

tx_lut <- full_join(
  coll_tx |>
    mutate(gene_id = str_match(unlist(Parent), "^GeneID:(WBGene\\d{8})$")[,2]) |>
    dplyr::select(gene_id,
                  tx_id = ID) |>
    distinct(),
  coll_exon |>
    mutate(tx_name = stringr::str_match(Name,"^([A-Za-z0-9.+_]+)\\.e[0-9]+$") %>%
             (function(x) x[,2])(),
           tx_id = unlist(Parent)) |>
    dplyr::select(tx_name, tx_id) |>
    distinct(),
  by = "tx_id")


write_csv(tx_lut,
          "data/intermediates_collapse/210910_tx_gene_LUT.csv")

