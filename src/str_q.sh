#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=str_q
#SBATCH -c 18
#SBATCH --mem=30G
#SBATCH --time=5-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# Principle: receive as argument a file with the sample names, and the name of the alignment pipeline
# *
# Pipeline: StringTie2 on results of bsn5 (alignment), using collapsed gtf,
# version 1, withOUT novel isoform discovery and GTF combination
# bsn5: bbduk, star, NuDup; then stc_dcq: stringtie on collapsed for novel, reference-based quantify
# Important: combine only same neuron
pipeline_version="str_q7"


######## Usage ######
#
# Call with:
#			sbatch str_q.sh
#
# Arguments:
# 	None.
#
# Note: if you modify this script, change the version number and create a new final directory (and if relevant re-run it on existing samples)
#
# Organization: takes the aligned files generated by another script, along with Wormbase GTF,
# then performs sample-wise quantification and saves results in output directory
# then calls prepDE.py and correct_gtf.gawk to generate the summary files (in output)
#
# 
# Inputs:
#     * alignments bam (e.g. bsn9)
#     * GTF (e.g. collapsed GTF or Wormbase GTF)
# 
# 
# Output:
# 			 * sample-wise quantifications (subdir quantifications)
# 			 * summary transcript and gene tables (subdir summaries)
# 
# 
#####################


echo "Starting $pipeline_version on $(date)"
echo



#-------------------------------------------------------------------------------------------------------------
# Create workspace paths
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references"

WSversion="WS281"

ref_gtf=$ref_dir/${WSversion}/"c_elegans.PRJNA13758."${WSversion}".canonical_geneset.gtf"

# use the merged alignments in scratch60 (technical replicates are merged)
alig_dir="/home/aw853/scratch60/2021-11-08_alignments"

str2_out="intermediates/220128_str_q_outs"


## Check inputs

if [ ! -d $alig_dir ]
then
echo "Error: bam directory does not exist: $alig_dir"
exit 1
fi


if [ ! -d $str2_out ]
then
echo "Error: destination directory does not exist. If this is a new pipeline version, you need to create it yourself: $str2_out"
exit 1
fi

mkdir -p $str2_out/neurons
mkdir -p $str2_out/quantifications
mkdir -p $str2_out/summaries


echo "--------------------------------------------------------------"
## Read sample list and remove trailing extension

echo "Reading samples from $alig_dir"
mapfile -t sampleList < <(ls $alig_dir/*.bam | xargs basename -a -s .bam)
mapfile -t neuronList < <(ls $alig_dir/*.bam | xargs basename -a -s .bam | cut -f1 -d"r")

if [ ${#sampleList[@]} -lt 1 ]
  then
  echo "Error: failed to find samples."
  exit 1
fi
  
if [ ${#neuronList[*]} -ne ${#sampleList[*]} ]
then
  echo "Error: number of samples and number of neurons do not match."
  exit 1
fi


# make lists of paths
for((i=0; i<${#sampleList[@]}; i++))
do
  samplePath[i]=$alig_dir/${sampleList[i]}".bam"
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


# list neurons
mapfile -t neurons < <(printf '%s\n' "${neuronList[@]}" | uniq)
nb_neurons=${#neurons[@]}

echo " Will treat ${#sampleList[@]]} samples from $nb_neurons neurons"









echo "------------------------------------------------------------------------------------------"
echo "--------------------------------    StringTie2 quantify    -------------------------------"
echo "------------------------------------------------------------------------------------------"
echo

for((i=0; i<$nb_samples; i++))
do
  sample=${sampleList[i]}
  
  echo "###################      Processing sample: $sample      ###################"
  echo
  
  mkdir -p $str2_out/quantifications/$sample
  stringtie2 -eB \
             -p $SLURM_CPUS_PER_TASK \
             -G $ref_gtf \
             -o $str2_out/quantifications/$sample/$sample.gtf \
             ${samplePath[i]}
done
echo

echo "-----------     Making tables     ------------"
echo
echo " > command: ~/.utilities/prepDE.py -i $str2_out"
cd $str2_out/summaries
~/.utilities/prepDE.py -i ../quantifications




echo
echo
echo "----------------------------------------------------------------------------------------------------------------------------------------------------------"
echo "Done: $(date)"