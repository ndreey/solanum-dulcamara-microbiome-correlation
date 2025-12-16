#!/bin/bash

#SBATCH --job-name PopStruct
#SBATCH -A naiss2025-22-494
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=36GB
#SBATCH -t 5:30:00
#SBATCH --output=slurm-logs/PopStruct/SLURM-%j.out
#SBATCH --error=slurm-logs/PopStruct/SLURM-%j.err
#SBATCH --mail-user=andbou95@gmail.com
#SBATCH --mail-type=ALL


# Start time and date
echo "$(date) [START]     Starting script execution"

# Load in modules
ml load bioinfo-tools
ml load bcftools/1.20
ml load plink/2.00a5.14
ml load ADMIXTURE/1.3.0


# Arguments
CPU=12
VCF=$1
ID=$(basename $VCF .vcf.gz)
WINDOW=$2
STEP=$3
R=$4
PREFIX=$ID-$WINDOW-$STEP-$R
LD_OUT=01-pruneLD/$PREFIX
PCA_OUT=02-PCA/$PREFIX
ADMIX_OUT=03-ADMIXTURE/$PREFIX
WORK=$(pwd)

mkdir -p $LD_OUT $PCA_OUT $ADMIX_OUT


echo "$(date) [INFO]     LD pruning"
plink2 \
    --vcf $VCF \
    --threads $CPU \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids "@:#" \
    --indep-pairwise $WINDOW $STEP $R \
    --out $LD_OUT/$PREFIX

echo "$(date) [EXEC]     Run PCA"
plink2 \
    --vcf $VCF \
    --threads $CPU \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids "@:#" \
    --extract $LD_OUT/$PREFIX.prune.in \
    --make-bed \
    --pca \
    --out $PCA_OUT/$PREFIX

echo "$(date) [Info]     Setting up for ADMIXTURE"
# ADMIXTURE does not accept chromosome names that are not human chromosomes.
# Change name of original
mv $PCA_OUT/$PREFIX.bim $PCA_OUT/$PREFIX.human-chr.bim

# Create a new one with first column set to 0
cat $PCA_OUT/$PREFIX.human-chr.bim | \
    awk '{$1="0"; OFS="\t"; print $0}' > $PCA_OUT/$PREFIX.bim

# Move to admix folder as $PREFIX.$K.P and $PREFIX.$K.Q are outputed
# to cwd.
cd $ADMIX_OUT

for K in {1..9}; do
    echo "$(date) [EXEC]     Running admixture for $K"
    admixture --cv -j8 ../../$PCA_OUT/$PREFIX.bed $K > $PREFIX-log$K.out
done

# Move back to working directory
cd $WORK

# Load in the r-arena mamba environment
source /cfs/klemming/home/a/andbou/.bashrc
mamba activate r-arena

# PCA specific arguments
EIGEN_VAL=$PCA_OUT/*.eigenval
EIGEN_VEC=$PCA_OUT/*.eigenvec

echo "$(date) [Info]     Plotting PCA"
Rscript scripts/r-scripts/plotPCA.R \
    $EIGEN_VAL \
    $EIGEN_VEC \
    $PCA_OUT

# Admix specific arguments
cat $PCA_OUT/*.fam | cut -f 1 > $PCA_OUT/samples.txt

ADMIX_PREFIX=$ADMIX_OUT/$PREFIX
ADMIX_SAMPLES=$PCA_OUT/samples.txt

echo "$(date) [Info]     Plotting ADMIXTURE"
for K in {2..9}; do
    echo "$(date) [Info]     K:$K"
    Rscript scripts/r-scripts/plotADMIXTURE.R \
        $ADMIX_PREFIX \
        $ADMIX_SAMPLES \
        $ADMIX_OUT \
        $K
done
echo "$(date) [COMPLETE]     Done!"