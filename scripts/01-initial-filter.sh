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
VCF=$1
PREFIX=$2
VCF_ANNO=vcfs/$PREFIX.raw.anno.vcf.gz
VCF_FILTERED=vcfs/$PREFIX.flagged.vcf.gz
VCF_OUT=vcfs/$PREFIX.filt.vcf.gz
CHR_RENAME=doc/chr_map_replace-invariants.txt

echo "$(date) [INFO]        Renaming chromosomes and annotating VCF with calculated tags"
# Get the chromosome names and create a file with <old_chr>white_space<new_name>
bcftools query -f '%CHROM\n' $VCF | sort -u | awk '{print $1,"rad_"$1}' > $CHR_RENAME

# Rename the chromosomes and fill out the tags.
bcftools annotate --threads $CPU --rename-chr $CHR_RENAME $VCF |
    bcftools +fill-tags --threads $CPU -Oz -o $VCF_ANNO -- -t AC,AN,MAF,HWE,ExcHet,F_MISSING
tabix $VCF_ANNO

echo "$(date) [INFO]     Apply individual genotype filtering to missing if below thresholds"
# Lets do further quality filtering
# Replace individual GT with missing if depth <3 or quality is < 20
bcftools filter --threads $CPU -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o $VCF_FILTERED $VCF_ANNO

echo "$(date) [INFO]        Index final VCF!"
tabix $VCF_FILTERED

echo "$(date) [INFO]     Filtering on site missingngess and MAF"
# Exclude records with missingness above 20% and MAF below 0.03.
bcftools filter --threads $CPU -e 'F_MISSING > 0.2 || MAF <= 0.03' $VCF_FILTERED \
| bcftools view --threads $CPU -f PASS -O z -o $VCF_OUT 

echo "$(date) [INFO]     Index vcf"
tabix $VCF_OUT

echo "$(date) [INFO]     Filtering complete!"
echo "$(date) [FINISH]     Calculate stats and plot!"

rm $VCF_FILTERED $VCF_FILTERED.tbi

# Calculate stats for created files
for f in $VCF_ANNO $VCF_OUT; do
    sbatch scripts/getSNPstats.sh $f
done 

echo "$(date) [FINISH]     All done!"