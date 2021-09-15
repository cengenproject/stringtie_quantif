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

# Find single-isoform genes
isof_per_gene <- cnts_long |> 
  group_by(gene_id) |>
  summarize(nb_isof = n_distinct(tx_name))
list_single_isof_genes <- isof_per_gene |>
  filter(nb_isof < 2) |>
  pull(gene_id)

# Find usage of major/minor isoforms
cnts_sample_prop <- cnts_long |>
  filter(! gene_id %in% list_single_isof_genes) |>
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



## Visualization ----

cnts_long <- readRDS("data/intermediates_visualization/210910_cnts_long.rds")




my_gene <- "ric-4"


cnts_long |>
  filter(gene_id == s2i(my_gene,gids)) |>
  group_by(sample_id) |>
  mutate(prop_sample = value/sum(value)) |>
  group_by(tx_name, neuron) |>
  summarize(prop = mean(prop_sample),
         sd = sd(prop_sample),
         .groups = "drop") |>
  ggplot(aes(x=neuron,y=prop, color = tx_name)) +
  geom_point(position=position_dodge(.2), pch=15,size=2) +
  geom_errorbar(aes(ymin=prop-sd,ymax=prop+sd),width=0,position=position_dodge(.2)) +
  theme_classic() +
  scale_y_continuous(limits=c(0,1), oob = scales::squish) +
  coord_flip() +
  labs(x=NULL,y=NULL,title=my_gene) +
  ggsci::scale_color_d3("category20")


cnts_long |>
  filter(gene_id == s2i(my_gene,gids)) |>
  group_by(tx_name, neuron) |>
  summarize(mean_cnt = mean(value),
            sd_cnt = sd(value),
            .groups = "drop") |>
  group_by(neuron) |>
  mutate(prop = mean_cnt/sum(mean_cnt),
         sd = sd_cnt/sum(mean_cnt)) |>
  ggplot(aes(x=neuron,y=prop, fill = tx_name)) +
  geom_col(position = position_stack()) +
  # geom_errorbar(aes(ymin=prop,ymax=prop+sd), position = position_stack()) +
  theme_classic() +
  coord_flip() +
  labs(x=NULL,y=NULL,title=my_gene) +
  # hues::scale_fill_iwanthue() +
  # ggsci::scale_fill_npg() +
  guides(fill=guide_legend(title="Transcript")) +
  scale_y_continuous(labels = scales::percent)


# Ordered

cnts_long |>
  filter(gene_id == s2i(my_gene,gids)) |>
  group_by(tx_name, neuron) |>
  summarize(mean_cnt = mean(value),
            sd_cnt = sd(value),
            .groups = "drop") |>
  group_by(neuron) |>
  mutate(prop = mean_cnt/sum(mean_cnt),
         sd = sd_cnt/sum(mean_cnt)) |>
  ungroup() |>
  mutate(neuron = fct_reorder(neuron, prop, min)) |>
  ggplot(aes(x=neuron,y=prop, fill = tx_name)) +
  geom_col(position = position_stack()) +
  # geom_errorbar(aes(ymin=prop,ymax=prop+sd), position = position_stack()) +
  theme_classic() +
  coord_flip() +
  labs(x=NULL,y=NULL,title=my_gene) +
  # hues::scale_fill_iwanthue() +
  # ggsci::scale_fill_npg() +
  guides(fill=guide_legend(title="Transcript")) +
  scale_y_continuous(labels = scales::percent)





ind_graph <- cnts_long |>
  filter(gene_id == s2i(my_gene,gids)) |>
  ggplot(aes(x=sample_id,y=value, fill = tx_name)) +
  geom_col(position = position_stack()) +
  theme_classic() +
  coord_flip() +
  # scale_y_continuous(trans = "log1p") +
  labs(x=NULL,y=NULL,title=my_gene)

plotly::ggplotly(ind_graph)


# By Amount, not Proportion
cnts_long |>
  filter(gene_id == s2i(my_gene,gids)) |>
  group_by(sample_id) |>
  group_by(tx_name, neuron) |>
  summarize(mean = mean(value),
            sd = sd(value),
            .groups = "drop") |>
  ggplot(aes(x=neuron,y=mean, color = tx_name)) +
  geom_point(position=position_dodge(.2), pch=15,size=2) +
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0,position=position_dodge(.2)) +
  theme_classic() +
  scale_y_log10() +
  # scale_y_continuous(limits=c(0,1), oob = scales::squish) +
  coord_flip() +
  labs(x=NULL,y=NULL,title=my_gene) +
  ggsci::scale_color_d3("category20")




# By single transcript

# boxplot
cnts_long |>
  filter(gene_id == s2i(my_gene,gids),
         tx_name == "B0350.2a.1") |>
  group_by(neuron) |>
  mutate(med_neur = median(value)) |>
  ungroup() |>
  mutate(neuron = fct_reorder(neuron, med_neur, max)) |>
  ggplot(aes(x = neuron, y = value)) +
  geom_boxplot(fill="grey90", color="grey60") +
  geom_point() +
  theme_classic() +
  # scale_y_continuous(trans = "log1p") +
  coord_flip()

#barplot of median
cnts_long |>
  filter(gene_id == s2i(my_gene,gids),
         tx_name == "C18H7.2a.1+b.1") |>
  group_by(neuron) |>
  mutate(med_neur = median(value)) |>
  ungroup() |>
  mutate(neuron = fct_reorder(neuron, med_neur, median)) |>
  ggplot(aes(x = neuron, y = value)) +
  # geom_boxplot(fill="grey90", color="grey60") +
  # geom_point() +
  geom_col(aes(y = med_neur)) +
  theme_classic() +
  # scale_y_continuous(trans = "log1p") +
  coord_flip()





xx <- c("ADA","ADE","ADF","ADL","AFD","AIA","AIB","AIM","AIN","AIY","AIZ","ALA","ALM","ALN","AQR","ASE","ASG",
        "ASH","ASI","ASJ","ASK","AUA","AVA","AVB","AVD","AVE","AVF","AVG","AVH","AVJ","AVK","AVL","AVM","AWA",
        "AWB","AWC","BAG","BDU","CAN","CEP","DVA","DVB","DVC","FLP","HSN","I1","I2","I3","I4","I5","I6","IL1",
        "IL1","IL2","IL2","LUA","M1","M2","M3","M4","M5","MC","MI","NSM","OLL","OLQ","PDA","PDB",
        "PDE","PHA","PHB","PHC","PLM","PLN","PQR","PVC","PVD","PVM","PVN","PVP","PVQ","PVR","PVT","PVW",
        "RIA","RIB","RIC","RID","RIF","RIG","RIH","RIM","RIP","RIR","RIS","RIV","RMD","RMD","RME",
        "RMF","RMG","RMH","SAA","SAB","SDQ","SIA","SIB","SMB","SMD","URA","URB","URX","URY","DD","VD","DA",
        "VA","DB","VB","AS","VC")

clipr::write_clip(xx %in% samples_table$neuron)

xx <- xx[xx %in% samples_table$neuron]


cnts_long |>
  filter(gene_id == s2i(my_gene,gids),
         tx_name == "C18H7.2b.2") |>
  group_by(neuron) |>
  summarize(mean_neur = mean(value)) |>
  ungroup() |>
  pull(mean_neur) |>
  round() |>
  clipr::write_clip()

cnts_long |>
  filter(gene_id == s2i(my_gene,gids),
         tx_name == "C18H7.2b.2") |>
  group_by(neuron) |>
  summarize(mean_neur = mean(value)) |>
  ungroup() %>%
  {set_names(round(.$mean_neur), .$neuron)}

