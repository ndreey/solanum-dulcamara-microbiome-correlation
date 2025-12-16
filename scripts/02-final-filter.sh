#!/bin/bash

#SBATCH --job-name sample-filter
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH -t 01:00:00
#SBATCH --output=slurm-logs/wrangle/SLURM-%j.out
#SBATCH --error=slurm-logs/wrangle/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL

# Start time and date
echo "$(date) [START]     Starting sample filtering script execution"

ml load bcftools/1.20

# Arguments
CPU=8
VCF=$1
FMISS_THRESH=0.3

# Define output files
PREFIX=$(basename $VCF .vcf.gz)
STATS_DIR=stats/$PREFIX
VCF_CLEAN=vcfs/${PREFIX}.clean.vcf.gz
VCF_CLEAN_ANNO=vcfs/${PREFIX}.clean.anno.vcf.gz
VCF_FLAG=vcfs/${PREFIX}.clean.flag.vcf.gz
VCF_OUT=vcfs/${PREFIX}.final.vcf.gz
SAMPLES_TO_REMOVE=stats/samples_to_remove.txt

# Find corresponding .imiss file
IMISS_FILE=$STATS_DIR/$PREFIX.imiss

echo "$(date) [INFO]     Checking for missingness file: $IMISS_FILE"

echo "$(date) [INFO]     Identifying samples with F_MISS > $FMISS_THRESH"
# Find samples with high missingness (skip header, column 5 is F_MISS)
HIGH_MISS_SAMPLES=$(awk -v thresh=$FMISS_THRESH 'NR>1 && $5 > thresh {print $1}' $IMISS_FILE)

if [[ -z "$HIGH_MISS_SAMPLES" ]]; then
    echo "$(date) [INFO]     No samples exceed F_MISS threshold of $FMISS_THRESH"
    echo "$(date) [INFO]     Copying original VCF as clean VCF"
    cp $VCF $VCF_OUT
    tabix -f $VCF_OUT
    exit 1
else
    # Count and display samples to remove
    SAMPLE_COUNT=$(echo "$HIGH_MISS_SAMPLES" | wc -w)
    echo "$(date) [INFO]     Found $SAMPLE_COUNT samples exceeding F_MISS threshold:"
    echo "$HIGH_MISS_SAMPLES" | tr ' ' '\n'
    
    # Create comma-separated list for bcftools
    #SAMPLES_TO_REMOVE=$(echo "$HIGH_MISS_SAMPLES" | tr ' ' ',')
    
    # Create a .txt file with sample names that are to be removed
    echo "$HIGH_MISS_SAMPLES" | tr ' ' '\n' > $SAMPLES_TO_REMOVE

    echo "$(date) [INFO]     Removing high-missingness samples from VCF"
    bcftools view -S ^${SAMPLES_TO_REMOVE} $VCF --threads $CPU -O z -o $VCF_CLEAN
    tabix $VCF_CLEAN
fi

echo "$(date) [INFO]     Re-annotating cleaned VCF with calculated tags"
bcftools +fill-tags $VCF_CLEAN --threads $CPU -Oz -o $VCF_CLEAN_ANNO -- -t AC,AN,MAF,HWE,ExcHet,F_MISSING
tabix $VCF_CLEAN_ANNO

echo "$(date) [INFO]     Apply individual genotype filtering to missing if below thresholds"
# Replace individual GT with missing if depth <3 or quality is < 20
bcftools filter --threads $CPU -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o $VCF_FLAG $VCF_CLEAN_ANNO

echo "$(date) [INFO]     Index final filtered VCF"
tabix $VCF_FLAG

echo "$(date) [INFO]     Filtering on site missingngess and MAF"
# Exclude records with missingness above 30% and MAF below 0.03.
bcftools filter --threads $CPU -e 'F_MISSING > 0.2 || MAF <= 0.03' $VCF_FLAG \
| bcftools view --threads $CPU -f PASS -O z -o $VCF_OUT 

echo "$(date) [INFO]     Index vcf"
tabix $VCF_OUT

echo "$(date) [INFO]     Filtering complete!"
echo "$(date) [FINISH]    Calculate stats and plot!"

rm $VCF_FLAG $VCF_CLEAN $VCF_FLAG.tbi $VCF_CLEAN.tbi

# Calculate stats for all VCF versions
for f in $VCF_CLEAN_ANNO $VCF_OUT; do
    if [[ -f "$f" ]]; then
        sbatch scripts/getSNPstats.sh $f
    fi
done

# Also run bcftools stats
sbatch scripts/getBCFstats.sh

echo "$(date) [FINISH]    All done!"