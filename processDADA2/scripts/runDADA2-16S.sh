#!/bin/bash

#SBATCH --job-name dada2
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=126GB
#SBATCH -t 08:30:00
#SBATCH --output=slurm-logs/dada2/SLURM-dada2-%j.out
#SBATCH --error=slurm-logs/dada2/SLURM-dada2-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in the mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate slu_proj_R


echo "$(date) [Info]     Run the DADA2 16S processing"
Rscript scripts/r-scripts/procesdada2.R


echo "$(date) [DONE]       Complete!"