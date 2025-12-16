#!/bin/bash

#SBATCH --job-name piawka
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=36GB
#SBATCH -t 05:30:00
#SBATCH --output=slurm-logs/piawka/SLURM-%j.out
#SBATCH --error=slurm-logs/piawka/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Load in the r-arena mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate vcf

# Arguments
CPU=8
JOBID=${SLURM_ARRAY_TASK_ID}
OUT_DIR=04-piawka
mkdir -p $OUT_DIR

VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
POP_MAP=doc/pop_map.tsv
POP_TMP=$(mktemp $OUT_DIR/$PREFIX-XXXXXX.samples)
POPS=$(mktemp $OUT_DIR/$PREFIX-XXXXXX.groups)
FMISS=0.3

echo "$(date) [INFO]        Mapping sample names to population"
# Get samples in VCF
bcftools query -l $VCF > $POP_TMP

# Make pop file
cat $POP_MAP | grep -f $POP_TMP > $POPS

echo "$(date) [INFO]        Running piawka"
# -M, --miss maximum share of missing genotypes at a site per group. 
# Higher-value sites are skipped. Should be a number between 0 and 1; default 1
piawka --vcf $VCF --groups $POPS --fst --jobs $CPU --miss $FMISS > $OUT_DIR/$PREFIX.tsv

echo "$(date) [FINISH]        Script all done!"