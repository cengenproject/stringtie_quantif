#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=run_R
#SBATCH -c 15
#SBATCH --mem=70G
#SBATCH --time=18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script $1 ---------------------"

module load R
R --slave -f R/collapse_tx.R

