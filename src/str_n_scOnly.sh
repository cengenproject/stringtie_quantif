#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=stringtie_novel_lr
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH --time=3-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# Principle: Rebuild txome on mix short-reads data.
# Same script as str_sc_n.sh, but not using long-reads: what is their contribution?
# *
# Pipeline:
#    * bulk short-reads aligned by STAR, all grouped in single BAM file (e.g. bsn9)
#       Note: actually reusing "subsampled_merged_short_reads.bam" generated on 2022-03-22
#    * NOT USING: 10x + IsoSeq long-reads, aligned by minimap2
#    * WS annotation
# bsn9: bbduk, star, NuDup; then str_n: stringtie on reference, novel discovery
pipeline_version="str_n"


######## Usage ######
#
# Call with:
#			sbatch str_n.sh
#
# Arguments:
# 	None.
#
# Note: if you modify this script, change the version number and create a new final directory (and if relevant re-run it on existing samples)
#
# 
# 
#####################



echo "Starting $(date)."


#-------------------------------------------------------------------------------------------------------------
# --- workspace paths
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references"

WSversion="WS281"

ref_gtf=$ref_dir/$WSversion/c_elegans.PRJNA13758.$WSversion.canonical_geneset.gtf


# --- outputs
str2_int="/gpfs/ycga/scratch60/ysm/hammarlund/aw853/2022-03-22_str_sc_n"


str2_out="intermediates/2022-05-13_str_n"



## Check inputs


if [ ! -d "$str2_out" ]
then
echo "Note: destination directory does not exist. Creating: $str2_out"
mkdir -p $str2_out
fi

mkdir -p $str2_out/quantifications
mkdir -p $str2_out/summaries


echo "--------------------------------------------------------------"
## Read sample list and remove trailing extension








module load SAMtools







echo "------------------------------  StringTie2 discovery   -----------------------------"
echo


stringtie2 -p $SLURM_CPUS_PER_TASK \
            -G $ref_gtf \
            -o $str2_out/220513_merged_L.gtf \
            -L $str2_int/lr_sorted.bam







echo
echo
echo "----------------------------------------------------------------------------------------------------------------------------------------------------------"
echo "Done $(date)."
echo


