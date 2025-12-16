#!/bin/bash

#SBATCH --job-name R_plot
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=56GB
#SBATCH -t 2:30:00
#SBATCH --output=slurm-logs/R/SLURM-%j.out
#SBATCH --error=slurm-logs/R/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in the r-arena mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate r-arena

VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
DIR=stats/$PREFIX
STATS_TSV=$DIR/raw_stats-$PREFIX.tsv
OUT_DIR=$DIR/plots

FRQ=$DIR/$PREFIX.frq
idepth=$DIR/$PREFIX.idepth
imiss=$DIR/$PREFIX.imiss
lmiss=$DIR/$PREFIX.lmiss

echo "$(date) [Info]     Plotting genome wide stats"
Rscript scripts/r-scripts/plotRawStats.R \
    $STATS_TSV \
    $OUT_DIR

echo "$(date) [Info]     Plotting vcftools stats"
Rscript scripts/r-scripts/plotVCFTOOLSstats.R \
    $FRQ $idepth $imiss $lmiss $OUT_DIR

echo "$(date) [DONE]       Complete!"