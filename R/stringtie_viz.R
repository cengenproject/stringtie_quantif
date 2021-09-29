## Initializations ----

suppressPackageStartupMessages({
  library(tidyverse)
})

library(wbData)
gids <- wb_load_gene_ids(277)






## Visualization ----

cnts_long <- readRDS("data/intermediates_visualization/210910_cnts_long.rds")




my_gene <- "nlg-1"


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
  # scale_y_log10() +
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

