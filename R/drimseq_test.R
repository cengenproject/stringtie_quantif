# Run after "drimseq_load_data"


## Initializations ----
library(tidyverse)
library(DRIMSeq)
library(wbData)
tx2g_tab <- wb_load_tx2gene(281)
gids <- wb_load_gene_ids(281)

fitdms <- qs::qread("data/intermediates_drimseq/2023-03-30_drimseq_fitdms.qs")


# counts for filtering
counts_per_gene <- counts(fitdms) |>
  pivot_longer(-c(gene_id, feature_id),
               names_to = "sample_id",
               values_to = "count") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}$")[,2]) |>
  group_by(neuron_id, gene_id) |>
  summarize(count = sum(count))


# Prepare contrasts ----

design_no_int <- model.matrix(~ 0 + neuron_id, data = DRIMSeq::samples(fitdms))

factors_in_design <- colnames(design_no_int)
contr_neurs <- map_chr(factors_in_design,
                        ~ paste0(.x,
                                 " - (",
                                 paste0(factors_in_design |> setdiff(.x),
                                        collapse = "+"),
                                 ") / ",
                                 length(factors_in_design) - 1)) |>
  set_names(str_match(factors_in_design, "^neuron_id([[:alnum:]]{2,4})$")[,2])



# Test head neurons ----

# ALA not sequenced

# ASER
tdms_ASER <- dmTest(fitdms, contrast = limma::makeContrasts(contr_neurs[["ASER"]], levels = design_no_int))


table(results(tdms_ASER)$adj_pvalue < 0.1, useNA = "ifany")

results(tdms_ASER) |> head()

results(tdms_ASER) |>
  filter(adj_pvalue < 0.1) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  arrange(desc(lr)) |>
  as_tibble()

genes_expr_in_ASER <- counts_per_gene |>
  filter(neuron_id == "ASER") |> 
  select(-neuron_id)


results(tdms_ASER) |>
  filter(adj_pvalue < 0.1) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  arrange(desc(lr)) |>
  as_tibble() |>
  left_join(genes_expr_in_ASER, by = "gene_id") |>
  clipr::write_clip()


# AVA

tdms_AVA <- dmTest(fitdms, contrast = limma::makeContrasts(contr_neurs[["AVA"]], levels = design_no_int))


table(results(tdms_AVA)$adj_pvalue < 0.1, useNA = "ifany")

results(tdms_AVA) |> head()

genes_expr_in_AVA <- counts_per_gene |>
  filter(neuron_id == "AVA") |> 
  select(-neuron_id)


results(tdms_AVA) |>
  filter(adj_pvalue < 0.1) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  arrange(desc(lr)) |>
  as_tibble() |>
  left_join(genes_expr_in_AVA, by = "gene_id") |>
  filter(count > 2000) |>
  clipr::write_clip()





## Old code below this point ----




myprop <- proportions(fitdms) %>%
  filter(gene_id == geneid)

predict_isof_switch <- function(dms, geneid, neur){
  # Takes the two transcripts with most variation and predict 
  # whether same isoform based on name
  myprop <- proportions(dms) %>%
    filter(gene_id == geneid)
  tibble(
    feature_id = myprop$feature_id,
    in_neur = rowMeans(myprop[,as.character(samples_table$col_name[samples_table$neuron == neur])]),
    others = rowMeans(myprop[,as.character(samples_table$col_name[samples_table$neuron != neur])])
  ) %>%
    mutate(diff = abs(in_neur - others)) %>%
    arrange(desc(diff)) %>%
    dplyr::slice(1:2) %>%
    pull(feature_id) %>%
    str_match("^[A-Z0-9]+\\.[0-9]{1,2}([a-z])\\.[0-9]+$") %>%
    magrittr::extract(3:4) %>%
    tibble::enframe() %>%
    tidyr::pivot_wider() %>%
    transmute(`1` == `2`) %>%
    as.character() %>%
    switch("TRUE" = "no isoform switch",
           "FALSE" = "isoform switch",
           "NA" = "single isoform"
    ) %>%
    return()
}
