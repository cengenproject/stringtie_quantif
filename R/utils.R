# Create an arbitrary character to be used as seqname for novel genes
create_seq_name <- function(x){
  paste0("X",toupper(as.hexmode(as.integer(x))), ".1")
}


# Create a character to be used as gene ID
create_gene_id <- function(x){
  paste0("STGene", str_pad(x, 6, pad = 0))
}
