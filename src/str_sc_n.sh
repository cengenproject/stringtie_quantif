#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=stringtie_novel
#SBATCH -c 15
#SBATCH --mem=60G
#SBATCH --time=3-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# Principle: Rebuild txome on mix short-reads and long-reads sc data.
# *
# Pipeline:
#    * bulk short-reads aligned by STAR, all grouped in single BAM file (e.g. bsn9)
#    * 10x + IsoSeq long-reads, aligned by minimap2
#    * WS annotation
# bsn9: bbduk, star, NuDup; then str_sc_n: stringtie on reference, with sc, novel
pipeline_version="str_sc_n"


######## Usage ######
#
# Call with:
#			sbatch str_sc_n.sh
#
# Arguments:
# 	None.
#
# Note: if you modify this script, change the version number and create a new final directory (and if relevant re-run it on existing samples)
#
# 
# 
#####################





#-------------------------------------------------------------------------------------------------------------
# --- workspace paths
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references"

WSversion="WS281"

ref_gtf=$ref_dir/$WSversion/c_elegans.PRJNA13758.$WSversion.canonical_geneset.gtf


# --- aligned reads
sr_alig_dir="/home/aw853/scratch60/2021-11-08_alignments"
lr_alig_sam="/SAY/standard/mh588-CC1100-MEDGEN/pacbio/2021-11-05/pb_squanti3/dedup.fasta.sam"


# --- outputs
str2_int="/gpfs/ycga/scratch60/ysm/hammarlund/aw853/"$(date +"%Y-%m-%d")"_"$pipeline_version
 str2_int="/gpfs/ycga/scratch60/ysm/hammarlund/aw853/2022-02-24_str_sc_n"

lr_sorted_bam=$str2_int/"2021-11-05_dedup.fasta.sorted.bam"

str2_out="intermediates/2022-02-24_str_sc_n"



## Check inputs

if [ ! -d "$sr_alig_dir" ]
then
  echo "Error: bam directory does not exist: $sr_alig_dir"
  exit 1
fi

if [ ! -f "$lr_alig_sam" ]
then
  echo "Error: long-reads bam file does not exist: $lr_alig_sam"
  exit 1
fi


if [ ! -d "$str2_out" ]
then
  echo "Note: destination directory does not exist. Creating: $str2_out"
  mkdir -p $str2_out
fi

mkdir -p $str2_out/quantifications
mkdir -p $str2_out/summaries


echo "--------------------------------------------------------------"
## Read sample list and remove trailing extension

echo "Reading samples from bsn9"
mapfile -t sampleList < <(ls $sr_alig_dir/*.bam | xargs basename -a -s .bam)


if [ ${#sampleList[@]} -lt 1 ]
  then
  echo "Error: failed to find samples."
  exit 1
fi




# make lists of paths
for((i=0; i<${#sampleList[@]}; i++))
do
  samplePath[i]=$sr_alig_dir/${sampleList[i]}".bam"
done


nb_samples=${#sampleList[@]}

if [ $nb_samples -ne ${#sampleList[*]} ]
then
  echo "Error while making list of samples."
  exit 1
fi

if [ ${#samplePath[@]} -ne $nb_samples ]
then
  echo "Error while building paths."
  exit 1
fi



echo " Will treat ${#sampleList[@]} samples."




mkdir -p $str2_int

module load SAMtools

echo "--------------------------  Merging short-reads samples   --------------------------"

samtools merge -@ $SLURM_CPUS_PER_TASK \
      -o $str2_int/merged_short_reads.bam \
      $(echo $sr_alig_dir/*.bam)
>>>>>>> 41080b0f344eef6ed4c4f24540babc3689ccad34

samtools view --bam \
  --subsample 0.1 \
  --subsample-seed 0 \
  -@ $SLURM_CPUS_PER_TASK \
  -o $str2_int/subsampled_merged_short_reads.bam \
  $str2_int/merged_short_reads.bam
  



echo "----------------------------  Sorting long-reads file   ----------------------------"

#samtools sort -@ $SLURM_CPUS_PER_TASK \
#        -o $lr_sorted_bam \
#        $lr_alig_sam

echo "------------------------------  StringTie2 discovery   -----------------------------"
echo


stringtie2 -p $SLURM_CPUS_PER_TASK \
            -G $ref_gtf \
            -o $str2_out/merged.gtf \
            --mix $str2_int/subsampled_merged_short_reads.bam $lr_sorted_bam








echo
echo "--------------------------------    StringTie2 quantify    -------------------------------"
echo

for((i=0; i<$nb_samples; i++))
do
  sample=${sampleList[i]}
  
  echo "###################      Quantify sample: $sample      ###################"
  echo
  
  mkdir -p $str2_out/quantifications/$sample
  stringtie2 -eB \
              -p $SLURM_CPUS_PER_TASK \
              -G $str2_out/merged_full.gtf \
              -o $str2_out/quantifications/$sample/$sample.gtf \
              ${samplePath[i]}
done

echo
echo "-----------     Making tables     ------------"
echo

echo " > command: ~/.utilities/prepDE.py -i $str2_out"
cur_dir=$(pwd)
cd $str2_out/summaries
~/.utilities/prepDE.py -i ../quantifications
cd $cur_dir



echo
echo
echo "----------------------------------------------------------------------------------------------------------------------------------------------------------"
echo "Done."
echo


