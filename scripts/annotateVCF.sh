#!/bin/bash

#SBATCH --job-name bcftools-annotate
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=36GB
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/annotate/SLURM-%j.out
#SBATCH --error=slurm-logs/annotate/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

ml load bcftools/1.20

# Arguments
CPU=8
VCF=$1
PREFIX=$2
VCF_ANNO=vcfs/$PREFIX.annotated.vcf.gz
CHR_RENAME=doc/chr_map_replace-invariants.txt

echo "$(date) [INFO]        Renaming chromosomes and annotating VCF with calculated tags"
# Get the chromosome names and create a file with <old_chr>white_space<new_name>
bcftools query -f '%CHROM\n' $VCF | sort -u | awk '{print $1,"rad_"$1}' > $CHR_RENAME

# Rename the chromosomes and fill out the tags.
bcftools annotate --threads $CPU --rename-chr $CHR_RENAME $VCF |
    bcftools +fill-tags --threads $CPU -Oz -o $VCF_ANNO -- -t AC,AN,MAF,HWE,ExcHet,F_MISSING
tabix $VCF_ANNO

echo "$(date) [Done]        Finished renaming chromosomes and annotating VCF with calculated tags"