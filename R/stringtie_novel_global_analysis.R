# Based on DRIMSeq preprocessing script
# First need to run STAR alignment (bsn5) + collapsedGTF + StringTie2 (stc_nq), with novel discovery



## Initializations ----

suppressPackageStartupMessages({
  library(tidyverse)
})

library(wbData)
gids <- wb_load_gene_ids(277)






## Prepare Data ----

# From StringTie2 (stc_nq)
path_data <- "data/210929_stc_nq_summaries/"

# get tx <-> gene table from "create_tx-gene_lookuptable.R": need to run separately
tx_table <- read_csv("data/210929_stc_nq_summaries/2021-10-14_tx_table.csv",
                     col_types = "ccccc")

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
  left_join(tx_table,
            by = "transcript_id")

cnts_long  <- cnts |>
  pivot_longer(-c("transcript_id","gene_id","ref_gene_id","ref_gene_seq","ref_gene_name","tx_name"),
               names_to = "sample_id") |>
  mutate(neuron = str_match(sample_id,"^([A-Zef0-9]{2,4})r\\d{1,3}$")[,2])



# saveRDS(cnts_long, "data/intermediates_visualization/211014_cnts_long.rds")
cnts_long <- readRDS("data/intermediates_visualization/211014_cnts_long.rds")



#~ Count reads to find best candidates ----

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


exons_novel <- rtracklayer::readGFF("data/210929_stc_nq_summaries/merged.gtf",
                                    filter = list(type="exon"),
                                    tags = c("transcript_id","gene_id"))




# Find multi-isoform genes
isof_per_gene <- cnts_long |> 
  select(gene_id, transcript_id) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_isof = n())

multi_isof_genes <- isof_per_gene |>
  filter(nb_isof > 1) |>
  pull(gene_id)

# Find single-exon genes
ex_per_gene <- exons_novel |>
  group_by(gene_id) |>
  summarize(nb_exons = n())

multi_exon_genes <- ex_per_gene |>
  filter(nb_exons > 1) |>
  pull(gene_id)

plot(eulerr::euler(list(multi_ex=multi_exon_genes, multi_isof=multi_isof_genes)))
table(multi_isof_genes %in% multi_exon_genes)
#> all multi-isoforms have several exons, there are multi-exonic genes with no alt splicing



# Find usage of major/minor isoforms
cnts_sample_prop <- cnts_long |>
  filter(gene_id %in% multi_isof_genes) |>
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
hist(cnts_major_minor_isoform$prop_second, breaks = 100)

ggplot(cnts_major_minor_isoform) +
  theme_classic() +
  geom_histogram(aes(x=prop_second), bins = 100,color="white") +
  stat_ecdf(aes(x = prop_second))

cnts_major_minor_isoform |>
  arrange(prop_second) |>
  mutate(cumu = cumsum(prop_second),
         cumu = cumu*sum(cnts_major_minor_isoform$prop_second == 1)/max(cumu)) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x=prop_second), bins = 400, color="white", fill = "grey")
  geom_line(aes(x=prop_second, y = cumu), size = 1) 
  scale_y_continuous(trans = "sqrt")
  
xx <- table(cnts_major_minor_isoform$prop_second)
yy <- xx[xx > 50]

cnts_major_minor_isoform |>
  filter(! prop_second %in% names(yy)) |>
  arrange(prop_second) |>
  mutate(cumu = cumsum(prop_second),
         cumu = cumu*sum(cnts_major_minor_isoform$prop_second > .97 &
                           cnts_major_minor_isoform$prop_second < 1)/max(cumu)) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x=prop_second), bins = 75, color="white", fill = "grey") +
  geom_line(aes(x=prop_second, y = cumu), size = 1)


cnts_major_minor_isoform |>
  filter(! prop_second %in% names(yy)) |>
  arrange(prop_second) |>
  mutate(x = seq(0, 1, length.out = n()),
         cumu = cumsum(prop_second > x),
         cumu = cumu/max(cumu),
         dc = c(0,diff(c(0,diff(cumu))))) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x=prop_second, y=after_stat(density)), bins = 100, color="white", fill = "grey") +
  stat_ecdf(aes(prop_second),size=1)+
  # geom_line(aes(x=prop_second, y = cumu), size = 1)
  geom_line(aes(x=x, y = dc/max(dc)), size = .5, color="darkred")



d <- cnts_major_minor_isoform |>
  filter(! prop_second %in% names(yy)) |>
  pull(prop_second) |>
  density()
summary(d)
d
plot(d)
plot(d$x,d$y,type="l")
plot(d$x, c(0,diff(d$y)), type = "l"); abline(h=0,col="grey",lty="dashed"); abline(v=.5,col="darkred",lty="dotted")
plot(d$x, c(0,0,diff(diff(d$y))), type = "l"); abline(h=0,col="grey",lty="dashed"); abline(v=.5,col="darkred",lty="dotted")

infl <- c(FALSE, diff(diff(diff(d$y))>0)!=0)
points(d$x[infl ], d$y[infl ], col="blue")

x <- cnts_major_minor_isoform |>
  filter(! prop_second %in% names(yy)) |>
  pull(prop_second)

m <- mean(x)
xx <- tibble(prop_second = x,
             diff = x-m,
             cumsum = cumsum(diff),
             dc = c(0,diff(cumsum)))

ggplot(xx) +
  theme_classic() +
  geom_smooth(aes(seq_along(cumsum), y=dc))


nrow(cnts_major_minor_isoform)
table(cnts_major_minor_isoform$prop_second >= .5)
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



# Illustrate

cnts_long|>
  filter(gene_id %in% c("MSTRG.10", "MSTRG.1003"),
         sample_id %in% c("AVMr226","AVMr227","AVMr228","AIYr65","AIYr66","AIYr67"),
         tx_name %in% c("Y74C9A.5.1","Y74C9A.4+Y74C9A.5+Y74C9A.9.s6",
                                "Y110A7A.13.2+1","Y110A7A.13+Y110A7A.14+Y110A7A.12.s1")) |>
  View()
cnts_sample_prop|>
  filter(gene_id %in% c("MSTRG.10", "MSTRG.1003"),
         sample_id %in% c("AVMr226","AVMr227","AVMr228","AIYr65","AIYr66","AIYr67"),
         tx_name %in% c("Y74C9A.5.1","Y74C9A.4+Y74C9A.5+Y74C9A.9.s6",
                        "Y110A7A.13.2+1","Y110A7A.13+Y110A7A.14+Y110A7A.12.s1")) |>
  View()

cnts_neuron_mean|>
  filter(gene_id %in% c("MSTRG.10", "MSTRG.1003")) |>
  View()

cnts_major_minor_isoform|>
  filter(gene_id %in% c("MSTRG.10", "MSTRG.1003")) |>
  View()
