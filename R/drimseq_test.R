# Run after "drimseq_load_data"


## Initializations ----
library(tidyverse)
library(DRIMSeq)
library(wbData)
tx2g_tab <- wb_load_tx2gene(281)


fitdms <- readRDS("data/intermediates_drimseq/2022-05-03_fitdms.rds")





# Prepare contrasts ----

design_no_int <- model.matrix(~ 0 + neuron_id, data = DRIMSeq::samples(fitdms))

factors_in_design <- colnames(design_no_int)
contr_neurs <- map_chr(factors_in_design,
                        ~ paste0(.x,
                                 " - (",
                                 paste0(factors_in_design |> setdiff(.x),
                                        collapse = "+"),
                                 ") / ",
                                 length(factors_in_design) - 1))


# Test ----

tdms_full <- dmTest(fitdms, contrast= my_contrasts, verbose = 1)


my_contrasts <- limma::makeContrasts(contr_neurs[1], levels=design_no_int)
tdms_AFD <- dmTest(fitdms, contrast= my_contrasts)
table(results(tdms_AFD)$adj_pvalue < 0.05)



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
