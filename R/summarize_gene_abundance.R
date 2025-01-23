# after running str_q which calls stringtie -A
# gather the gene_abundance.tab files and save a single gene_TPMs.tsv


library(tidyverse)


qdir <- "intermediates/250122_str_q_outs/quantifications"
summdir <- "intermediates/250122_str_q_outs/summaries"

samples_table <- tibble(sample_id = list.dirs(qdir,full.names=FALSE,recursive = FALSE),
                        gene_path = file.path(qdir, sample_id, "gene_abundance.tab"))




samples_table$content <- purrr::pmap(
  samples_table,
  \(sample_id,gene_path){
    read_tsv(gene_path,
             col_types = cols(
               `Gene ID` = col_character(),
               `Gene Name` = col_character(),
               Reference = col_character(),
               Strand = col_character(),
               Start = col_double(),
               End = col_double(),
               Coverage = col_double(),
               FPKM = col_double(),
               TPM = col_double()
             )) |>
      rename(gene_id = `Gene ID`,
             gene_name = `Gene Name`) |>
      add_column(sample_id)
  }
)

samples_table$gene_path <- NULL


# all samples have the same genes (as we used ref annotation)
stopifnot(length(unique(map_int(samples_table$content, nrow))) == 1L)
stopifnot(all(map_int(samples_table$content, ncol) == 10L))


genes_long <- samples_table$content |>
  bind_rows() |>
  select(sample_id, gene_id, gene_name, TPM)

write_tsv(genes_long,
          file.path(summdir, "250123_gene_TPMs.tsv"))









