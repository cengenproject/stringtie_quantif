# Based on DRIMSeq preprocessing script
# First need to run STAR alignment (bsn5) + collapsedGTF + StringTie2 (stc_dq), no novel discovery



## Initializations ----

suppressPackageStartupMessages({
  library(tidyverse)
})

library(wbData)
gids <- wb_load_gene_ids(277)






## Prepare Data ----

# From StringTie2 (stc_dcq5)
path_data <- "data/210909_str2_summaries/"

# get tx <-> gene table from "find_gene_for_tx5.R": need to run separately
tx_lut <- read_csv("data/intermediates_collapse/210910_tx_gene_LUT.csv",
                   col_types = "ccc")

# Prepare samples table
samples_table <- readr::read_csv("data/210820_full_sample_list.csv", col_types = "ccc") %>%
  dplyr::rename(sample_id = sample)


# Data from STAR alignment + collapsedGTF + StringTie2 (bsn5cst1_1 pipeline)
# Read sample counts
cnts <- readr::read_csv(file.path(path_data, "transcript_count_matrix.csv"),
                        col_types = paste0(c("c", rep("d", nrow(samples_table))),
                                           collapse=""))


# Check samples
str_match(colnames(cnts), "^([A-Z0-9]{1,4}|Ref)r[0-9]{1,4}$") %>%
  as_tibble(.name_repair = "universal") %>%
  dplyr::slice(-1) %>%
  dplyr::rename(sample_id = ...1,
                neuron = ...2) %>%
  arrange(sample_id) %>%
  identical(samples_table %>% select(-replicate)) %>%
  stopifnot()


# rename genes and transcripts
cnts <- cnts |>
  left_join(tx_lut,
            by = c(transcript_id = "tx_id"))

cnts_long  <- cnts |>
  pivot_longer(-c("transcript_id","gene_id","tx_name"),
               names_to = "sample_id") |>
  mutate(neuron = str_match(sample_id,"^([A-Zef0-9]{2,4})r\\d{1,3}$")[,2])

# saveRDS(cnts_long, "data/intermediates_visualization/210910_cnts_long.rds")


#~ Count reads ----

candidates <- read_tsv("data/intermediates_visualization/candidates.txt",
                       col_names = "gene_name", col_types = "c")


cmat <- as.matrix(cnts[,2:170])

left_join(candidates |> 
            mutate(gene_id = s2i(gene_name, gids)),
          tibble(gene_id   = cnts$gene_id,
                 sum_count = apply(cmat, 1, sum)/1e5) |> 
            group_by(gene_id) |>
            summarize(sum_count = sum(sum_count))) |>
  pull(sum_count) |>
  clipr::write_clip()


#~ Count genes with AS ----

# Keep protein-coding, multi-exon genes
gids |> 
  filter(biotype == "protein_coding_gene")

exons <- wb_load_exon_coords(277)

ex_per_gene <- exons |>
  group_by(gene_id) |>
  summarize(nb_exons = n(),
            biotype = first(gene_biotype))

multi_ex_pcd <- ex_per_gene |> 
  filter(nb_exons > 1,
         biotype == "protein_coding") |> 
  pull(gene_id)

# Find single-isoform genes
isof_per_gene <- cnts_long |> 
  group_by(gene_id) |>
  summarize(nb_isof = n_distinct(tx_name))
list_single_isof_genes <- isof_per_gene |>
  filter(nb_isof < 2) |>
  pull(gene_id)

# Find usage of major/minor isoforms
cnts_sample_prop <- cnts_long |>
  filter(! gene_id %in% list_single_isof_genes,
         gene_id %in% multi_ex_pcd) |>
  group_by(gene_id, sample_id) |>
  mutate(prop_sample = value/sum(value)) |>
  filter(! is.na(prop_sample))

cnts_neuron_mean <- cnts_sample_prop |>
  group_by(gene_id, tx_name, neuron) |>
  summarize(prop = mean(prop_sample),
            sd = sd(prop_sample),
            .groups = "drop")

cnts_major_minor_isoform <- cnts_neuron_mean |>
  group_by(gene_id, tx_name) |>
  summarize(max_prop = max(prop),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(prop_highest = max(max_prop),
            prop_second = sort(max_prop, decreasing = TRUE)[[2]])

hist(cnts_major_minor_isoform$prop_highest, breaks = 100)
nrow(cnts_major_minor_isoform)
table(cnts_major_minor_isoform$prop_second > .15)
# as per Wang (2008), count a gene AS if minor isof > 15%
# https://www.nature.com/articles/nature07509

# Just protein-coding genes

cnts_long |> select(gene_id) |> distinct() |> left_join(gids) |> pull(biotype) |> table()
#> 19,995 prot-coding genes in data
left_join(cnts_major_minor_isoform, gids, by = "gene_id") |>
  filter(biotype == "protein_coding_gene") |>
  pull(prop_second) %>%
  {table(. > .15)}
#> 6,175 of the prot-coding genes have minor isof > 15%



