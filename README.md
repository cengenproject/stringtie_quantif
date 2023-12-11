

Use `str_q.sh` to quantify transcript usage (requires 3 GB, 14h). Inputs: bams, gtf.

With bsn12, version str_q8.


# Contents of `str_q`
First runs Stringtie `-eB` and saves the results of individual samples in `intermediates/231208_str_q_outs/quantifications`.

Then calls `R/summarize_stringtie_q_output.R` which compiles matrices of coverage, FPKM, TPM, and saves them in `.../summaries`.

In particular, saves the TPMs also in a long format, as `.../t_exp.tsv`, for use with Shiny app. The SLURM log should give the scp command to transfer this file (see repo `isoforms_compare`).

Finally, calls `~/.utilities/prepDE.py3`, provided by Stringtie authors, to save `.../summaries/gene_count_matrix.csv` and `transcripts_count_matrix.csv` for other downstream uses (such as `R/drimseq_load_data` and `R/drimseq_test`). Note: 231208 had to change to the Python 3 version on McCleary.

## DRIMSeq

No definitive approach at this point, see `R/drimseq_load_data` to load and pre-filter, saving the object `intermediates/2023-03-30_drimseq_fitdms.qs`, that can be loaded in `R/drimseq_test` to perform DTU.



# Old approaches

Older versions used gtf augmented with novel transcripts (ref_gtf="intermediates/2022-03-22_str_sc_n/220322_novel_filt_sorted.gtf")


str_sc_n: novel isoforms using mix of short and long reads
str_n: novel isoforms only with short reads (to compare with `str_sc_n`)
str_n_scOnly: novel isoforms only with long reads (to compare with `str_sc_n`)


Older versions needed to prepare bam files with `prep_alignments.sh` (and save them in scratch60).
