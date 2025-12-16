#!/bin/bash

#SBATCH --job-name bcftools-wrangle
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=36GB
#SBATCH -t 00:30:00
#SBATCH --output=slurm-logs/wrangle/SLURM-%j.out
#SBATCH --error=slurm-logs/wrangle/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting script execution"

ml load bcftools/1.20

# Arguments
CPU=8
VCF_BIAL=vcfs/new.filt.final.vcf.gz
VCF_ALL=vcfs/populations.all.annotated.vcf.gz
VCF_INVAR=vcfs/final.invariants.vcf.gz
VCF_MERGED=vcfs/final.merged.vcf.gz
CHR_RENAME=doc/chr_map_replace-invariants.txt
SAMPLES_TO_REMOVE=stats/samples_to_remove.txt

echo "$(date) [INFO]        Renaming chromosomes and annotating VCF with calculated tags"
# Get the chromosome names and create a file with <old_chr>white_space<new_name>
bcftools query -f '%CHROM\n' $VCF_ALL | sort -u | awk '{print $1,"rad_"$1}' > $CHR_RENAME

# Rename the chromosomes and keep only invariants
bcftools annotate --threads $CPU --rename-chrs $CHR_RENAME $VCF_ALL | \
    bcftools view --threads $CPU -C 0 | \
    bcftools view --threads $CPU -S ^${SAMPLES_TO_REMOVE} -O z -o $VCF_INVAR
tabix $VCF_INVAR

echo "$(date) [INFO]     Merging biallelic SNPs with invariant sites"
# Merge the biallelic SNPs with invariant sites
bcftools concat --threads $CPU --allow-overlaps \
    $VCF_BIAL $VCF_INVAR | \
    bcftools sort -O z -o $VCF_MERGED
tabix $VCF_MERGED

# Also run bcftools stats
sbatch scripts/getBCFstats.sh

echo "$(date) [INFO]     Checking merged VCF statistics"
# Quick stats check
echo "Original biallelic SNPs:"
bcftools stats $VCF_BIAL | grep "number of records"

echo "Invariant sites:"
bcftools stats $VCF_INVAR | grep "number of records"

echo "Merged VCF:"
bcftools stats $VCF_MERGED | grep "number of records"

echo "$(date) [FINISH]   Script complete!"
