# Export the list of main transcripts for each gene-neuron, in a format to be used
# for the Deep Learning model.
# 
# Use the result of "DL_select_main_tx" which was run on the HPC to to the main computation
# and save an intermediate, load the intermediate, remove the genes that are not 
# expressed based on scRNA-Seq thresholded data, and export in a tsv.
# Ensure it matches the tx structures generated separately, see also the validity tests.


## Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

#~ load data ----

tx_neurons <- qs::qread("data/intermediates_for_DL/220919_tx_neuron.qs") |>
  ungroup()

gene_expression <- read.csv("data/thresholded_gene_expression/bsn9_subtracted_integrated_binarized_expression_withVDDD_FDR_0.1_092022.tsv",
                            sep = "\t") |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expression")



# Clean up: if not present in gene_expression, almost surely not expressed in any neuron
tx_neurons_filt <- tx_neurons |>
  filter(neuron_id %in% gene_expression$neuron_id,
         gene_id %in% gene_expression$gene_id)


tx_neurons_thresholded <- full_join(tx_neurons_filt,
                                    gene_expression,
                                    by = c("gene_id", "neuron_id")) |>
  filter(expression == 1L)

tx_neurons_thresholded |>
  select(gene_id,
         neuron_id,
         main_transcript = main_tx) |>
  write_tsv("data/intermediates_for_DL/220920_main_transcript_per_neuron.tsv")


