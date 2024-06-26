#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=merge_index_bam
#SBATCH -c 10
#SBATCH --mem=25G
#SBATCH --time=1-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# Transfer the bam files merging the technical replicates. Filter out multimappers. Index the results (bai)

echo "Transfer, merge, filter multimapper, index; starting $(date)"

module load SAMtools


alig_dir_orig="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_bams"
alig_dir="/home/aw853/scratch60/2022-03-18_alignments"

mkdir -p $alig_dir

mapfile -t sampleList < <(ls $alig_dir_orig/*.bam \
                          | xargs basename -a -s .bam \
                          | cut -f1 -d't'\
                          | uniq)

echo ${#sampleList[@]}" samples"

echo "Merging and filtering."

for sample in ${sampleList[@]}
do
  samtools merge -@ $SLURM_CPUS_PER_TASK -o - $(echo $alig_dir_orig/$sample"*.bam") \
  | samtools view -b --tag NH:1 \
  > $alig_dir/$sample".bam"
done

echo "BAM files merged and filtered. Indexing."

for sample in ${sampleList[@]}
do
	if [ -f "$alig_dir/$sample.bam.bai" ]; then
		echo "$sample already indexed."
	else 
		echo "Indexing $sample"
		samtools index -@ $SLURM_CPUS_PER_TASK $alig_dir/$sample".bam" $alig_dir/$sample".bam.bai"
	fi
done


echo
echo "All done: $(date)"
