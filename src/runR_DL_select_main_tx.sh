#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=run_R
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --time=2:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script 'DL_select_main_tx.R' ---------------------"

module load R
R --slave -f R/DL_select_main_tx.R

