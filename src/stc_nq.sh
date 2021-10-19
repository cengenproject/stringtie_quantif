#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=stringtie
#SBATCH -c 18
#SBATCH --mem=30G
#SBATCH --time=5-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# Principle: receive as argument a file with the sample names, and the name of the alignment pipeline
# *
# Pipeline: StringTie2 on results of bsn5 (alignment), using collapsed gtf, version 1, with novel isoform discovery and GTF combination
# bsn5: bbduk, star, NuDup; then stc_dcq: stringtie on collapsed for novel, combine, quantify
# Important: combine only same neuron
pipeline_version="stc_nq"


######## Usage ######
#
# Call with:
#			sbatch stc_nq.sh
#
# Arguments:
# 	None.
#
# Note: if you modify this script, change the version number and create a new final directory (and if relevant re-run it on existing samples)
#
# Organization: takes the aligned files generated by another script,
# performs first round of transcript discovery (StringTie) and saves sample-wise gtf in intermediary
# then merges for each neuron type and saves both mergelist (list of samples from that neuron) and merged gtf in output directory
# then merges all the neuron-wise gtf into a full gtf (stored in intermediary with list of neuron-wise gtf to merge)
# then performs samples-wise quantification using merged_full.gtf and saves results in output directory
# then calls prepDE.py to generate the summary files (in output)
#
# 
# To re-run with a few more samples, do sample-wise discovery for these new samples, merge with existing neuron-wise gtf
# then regenerate the full merged gtf with all neurons (modified or not), and redo the quantifications for all samples
# as well as the summaries
# 
# Inputs:
#     * alignments bam (e.g. bsn5)
#     * GTF (e.g. collapsed GTF)
# 
# In intermediary (on scratch60):
#           * sample-wise (discovery) gtf,
# 					* full merged gtf
# 
# In output:
#        * neuron-wise mergelists (subdir neurons)
# 			 * neuron-wise merged gtf (subdir neurons)
# 			 * mergelist for full merge (subdir quantifications)
# 			 * sample-wise quantifications (subdir quantifications)
# 			 * summary transcript and gene tables (subdir summaries)
# 			 * summary corrected gtf (with transcript names as in table) (subdir summaries)
# 
# 
#####################





#-------------------------------------------------------------------------------------------------------------
# Create workspace paths
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references"

WSversion="WS277"

ref_gtf="intermediates/210909_collapsed_annotation_${WSversion}.gff3"

# use the merged alignments in scratch60 (technical replicates are merged)
alig_dir="/home/aw853/scratch60/2021-08-18_alignments"

#str2_int="/gpfs/ycga/scratch60/ysm/hammarlund/aw853/"$(date +"%Y-%m-%d")"_"$pipeline_version
str2_int="/gpfs/ycga/scratch60/ysm/hammarlund/aw853/2021-10-19_stc_nq"
str2_out="intermediates/211019_str2_stc_nq_big_gap"


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

echo "Reading samples from bsn5"
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
echo "---------------------------------  StringTie2 discovery   -------------------------------"
echo "------------------------------------------------------------------------------------------"
echo
echo

mkdir -p $str2_int
for ((i=0; i<$nb_samples; i++))
do
  sample=${sampleList[i]}
  
  echo "######################     Processing sample: $sample      ######################"
  echo
  
  # First run for each sample without eB, to create the GTFs and discover novel transcripts
  
  # -f  minimum isoform abundance (default 0.01)
  # -j  min nb of jction-covering reads (def 1)
  # -c  minimum read coverage allowed anywhere in transcript (def 1)
  stringtie2 -p $SLURM_CPUS_PER_TASK \
    -G $ref_gtf \
    -o $str2_int/$sample.gtf \
    -f 0.05 \
    -j 10 \
    -c 1.5 \
    -g 100 \
    ${samplePath[i]}
done

echo "------------------------------------------------------------------------------------------"
echo "--------------------------------     StringTie2 merge     --------------------------------"
echo "------------------------------------------------------------------------------------------"
echo
echo

echo ">  Initialize mergelists"
touch $str2_out"/quantifications/mergelist_full.txt"

for neur in ${neurons[@]}
do
  touch $str2_out"/neurons/mergelist_"$neur".txt"
  echo $str2_out/neurons/merged_$neur.gtf >> $str2_out/quantifications/mergelist_full.txt
done

echo
echo ">  Fill mergelists"

for((i=0; i<$nb_samples; i++))
do
  echo $str2_int"/"${sampleList[i]}".gtf" >> $str2_out/neurons/mergelist_${neuronList[i]}.txt
done

echo
echo "Perform merging of GTFs"

for neur in ${neurons[@]}
do
  echo "######################  Merge neuron: $neur ######################"
  echo
  
  # -m min length of transcript (def 50)
  # -c min coverage (def 0)
  # -T min TPM (def 0)
  # -f min fraction of isoform (def 0.01)
  stringtie2 --merge \
    -G $ref_gtf \
    -o $str2_out/neurons/merged_$neur.gtf \
    -p $SLURM_CPUS_PER_TASK \
    -m 50 \
    -c 1 \
    -T 2 \
    -f 0.10 \
    -g 400 \
    $str2_out/neurons/mergelist_$neur.txt
done

echo
echo "Make master merge"


# -m min length of transcript (def 50)
# -c min coverage (def 0)
# -T min TPM (def 0)
# -f min fraction of isoform (def 0.01)
stringtie2 --merge \
  -G $ref_gtf \
  -o $str2_out/merged_full.gtf \
  -p $SLURM_CPUS_PER_TASK \
  -m 50 \
  -c 1 \
  -T 2 \
  -f 0.01 \
  -g 400 \
  $str2_out/quantifications/mergelist_full.txt




echo
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