

# Background


## Biology

Each neuron expresses a set of genes. A gene has a primary sequence that is processed by splicing to give individual transcripts. Each neuron type uses a subset of the possible transcripts.

During splicing, some parts of the primary sequence are removed (introns), others are kept (exons). Because of alternative splicing, an intron in one neuron can be an exon in another.

## Files

I am providing 3 sets of files:

* gene sequence
* transcripts structures
* transcript expression

## Identifiers

Each gene has a unique identifier in the format: `WBGene12345678`.

Each transcript has a unique identifier string. The format is similar to `A00A0.1a.1`, but there are variations on that format, do not make assumptions. While the exact string conveys biological information, this doesn't matter here.

<details>
  <summary>in case the details are of interest</summary>
  <p>In `Y6B3B.5a.2`, `Y6B3B` is a region of the genome, `Y6B3B.5` is the 5th gene in that region, `Y6B3B.5a` is the first protein isoform for that gene, `Y6B3B.5a.2` is the second transcript for that isoform.
  So `Y6B3B.5b.1` is another isoform of the same gene, i.e. the coding sequence is different, whereas `Y6B3B.5a.1` has the same coding sequence but different transcribed regions</p>
</details>

<br/>


Each transcript is separated into segments, which can be either an exon or intron. They are described by their start and end positions in the gene coordinates.


# Gene sequence

This is almost the same file as previously, after filtering to keep only relevant genes (see below). In FASTA format, each record (starting with `>`) corresponds to a gene.


# Transcripts structure

The segments are described by their start and end, in gene coordinates, inclusive. For example, if a transcript looks like:

```
  exon    intron    exon
ACGCCGCG ATATATTA CGCAGCAGT
|   |     |    |     |    |
1   5    10   15    20   25
```

this will become:

```
start  end    type
------------------
    1    8    exon
    9   16  intron
   17   25    exon
```


**Important note**: this is in gene coordinates, so if the same gene has several transcripts that do not start from the same position, there can be a gap at the beginning! (or same at the end)

For example, this situation can happen:

```
>gene_sequence
ATCTAGA ACGCCGCG ATATATTA CGCAGCAGT TGTCTTAAAC
|   |     |    |     |     |    |     |    | 
1   5    10   15    20    25   30    35   40

>transcript1
|   |      exon | intron |   exon     |    |
|   |   ACGCCGCG ATATATTA CGCAGCAGT   |    |
|   |     |    |     |     |    |     |    | 
1   5    10   15    20    25   30    35   40

>transcript2
 exon  |    intron       |  exon 
ATCTAGA ACGCCGCG ATATATTA CGCAGCAGT TGTCTTAAAC
|   |     |    |     |     |    |     |    | 
1   5    10   15    20    25   30    35   40
```

this will become:

```
transcript_id  start  end    type
---------------------------------
  transcript1      8   15    exon
  transcript1     16   23  intron
  transcript1     24   32    exon
  
  transcript2      1    7    exon
  transcript2      8   23  intron
  transcript2     24   42    exon
```


# Transcript expression

First, for each sample and each gene, StringTie was run, that establishes a splicing graph and uses a flow maximization algorithm to quantify each transcript expression.

Then, for each gene in each sample, I selected the highest transcript, and for each neuron type I selected the majority vote among samples, thus getting for each gene and neuron type, a single "main transcript".

The result looks like:

```
       gene_id   neuron_id   main_transcript
WBGene00000001         ADL           A01A1.1
WBGene00000001         PVD           A01A1.1
WBGene99999999         ADL        C07C5.3a.1
WBGene99999999         PVD        C07C5.3b.2
```

**Important note:** not every gene is expressed in every neuron type; if a gene is not expressed it is meaningless to look at its main transcript. So I filtered out the unexpressed genes based on previous scRNA-Seq data. So the result can look like:

```
       gene_id   neuron_id   main_transcript
WBGene00000001         ADL           A01A1.1
WBGene00000001         PVD           A01A1.1
WBGene00000002         ADL           B99B9.9
WBGene99999999         ADL        C07C5.3a.1
```
i.e. in that example the genes `WBGene00000002` and `WBGene99999999` are not expressed in neuron `PVD`, so there is no corresponding row.


# Datasets sizes

After thresholding, we have `10,712` genes, in `46` neuron types, so that `14,469` transcripts can be the main transcript in at least one neuron.

Thus, there are up to `10,712*46 = 492,752` neuron-gene combinations, but in practice only  `174,298` after removing the neuron-gene combinations that are not expressed.

There are `293,515` individual segments (exon or intron) across all genes.

Note: before thresholding, there were 42-46k genes, thus previous datasets would have that many genes (e.g. `42,075` genes in the FASTA file). The many genes removed are not very relevant for this task, most are not protein coding, others are not expressed in neurons etc...



# Datasets versions

StringTie quantification: `220322 str_q`.

scRNA-Seq: from Alec Barrett, `bsn9 subtracted integrated binarized expression withVDDD FDR 0.1 092022`, TPR=74%, FDR=10%.

select_main_tx: `220919`


<br/><br/><br/>

