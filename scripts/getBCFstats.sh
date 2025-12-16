#!/bin/bash

#SBATCH --job-name bcftools-stats
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12GB
#SBATCH -t 02:30:00
#SBATCH --output=slurm-logs/stats/SLURM-bcftools-%j.out
#SBATCH --error=slurm-logs/stats/SLURM-bcftools-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date) [START]     Starting script execution"

# Laad modules
ml load vcftools/0.1.16
ml load bcftools/1.20

BCF_OUT=stats/vcf_stats.tsv


#########################################################################
#                  Summarize BCFTOOLS Stats
#########################################################################
# Get quick stats for the files
echo -e "$(date) [INFO]        Summarising stats from each VCF file: $BCF_OUT"
echo -e "VCF\tsamples\trecords\tinvariants\tSNPs\tMNPs\tindels\tothers\tmulti-sites\tmulti-SNPs" > $BCF_OUT

for FILE in $(ls vcfs/*.vcf.gz); do

    echo -e "$(date) [INFO]        Summarising $FILE"
    echo -e "$(basename $FILE)\t$(bcftools stats $FILE | grep -v "# SN," | grep -A 9 "# SN" | grep -v "#" | cut -f 4 | paste -sd$'\t')" >> $BCF_OUT

done

# Preview the file
echo -e ">> BCFTOOLS STATS DESCRIPTION << \n"
echo "number of records   .. number of data rows in the VCF"
echo "number of no-ALTs   .. reference-only sites, ALT is either '.' or identical to REF"
echo "number of SNPs      .. number of rows with a SNP"
echo "number of MNPs      .. number of rows with a MNP, such as CC>TT"
echo "number of indels    .. number of rows with an indel"
echo "number of others    .. number of rows with other type, for example a symbolic allele or a complex substitution, such as ACT>TCGA"
echo "number of multiallelic sites     .. number of rows with multiple alternate alleles"
echo "number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs"
echo "Note that rows containing multiple types will be counted multiple times, in each counter. For example, a row with a SNP and an indel increments both the SNP and the indel counter."
echo ""

# Print out the stats in a clean table
echo -e "\n >> BCFTOOLS STATS <<"
column -t "$BCF_OUT"


echo -e "\n$(date) [FINISH]       Script Complete!\n"
