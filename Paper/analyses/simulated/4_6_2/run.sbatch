#!/bin/sh

#SBATCH --job-name=run
#SBATCH --output=run.out
#SBATCH --account=pi-joet1
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load R
Rscript bg_avg_prediction.R
Rscript ref_free_prediction.R
