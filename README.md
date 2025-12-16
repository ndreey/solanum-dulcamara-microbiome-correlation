# Solanum dulcamara Population Genomics and Microbiome Analysis

This project investigates the relationship between host plant population genetics and root-associated microbiome composition in *Solanum dulcamara* across natural populations. Using RADseq data and 16S/ITS amplicon sequencing, we test whether host genetic distance (FST) correlates with microbiome beta-diversity (Bray-Curtis dissimilarity) through Mantel correlation tests. The analysis pipeline includes SNP calling and filtering, population structure analysis (PCA, admixture), calculation of genetic differentiation statistics (FST, nucleotide diversity), and microbiome community profiling using DADA2. We compare bacterial (16S), fungal (ITS), and arbuscular mycorrhizal fungi (AMF) communities across taxonomic levels to determine which microbial groups show the strongest correlation with host genetic structure, providing insights into the eco-evolutionary dynamics of plant-microbiome associations.

## Workflow

### 1. Initial VCF Filtering
Renames chromosomes, annotates VCF with quality metrics (AC, AN, MAF, HWE, etc.), filters individual genotypes with low depth (<3) or quality (<20), and removes sites with >20% missingness or MAF ≤0.03. The arguments it takes is path/to/vcf and a string to be prefix

```bash
# Filter raw SNPs with quality thresholds
sbatch scripts/01-initial-filter.sh original_vcfs/pe_pop3.snps.vcf.gz new
```

### 2. Annotate VCF Files
Annotation is required before calculating statistics with bcftools/vcftools. The arguments it takes is path/to/vcf and a string to be prefix

```bash
# Annotate original files
sbatch scripts/annotateVCF.sh original_vcfs/populations.all.vcf.gz populations.all
sbatch scripts/annotateVCF.sh original_vcfs/pe_pop3.snps.vcf.gz pe_pop3

# Standardize filtered VCF to match continuous chromosome blocks
bcftools view original_vcfs/pe_pop3.snps.vcf.gz \
    -R original_vcfs/pe_pop3_filt.snps.vcf.gz \
    -O z -o original_vcfs/pe_pop3_filt.standardized.vcf.gz

# Create index
tabix original_vcfs/pe_pop3_filt.standardized.vcf.gz 

# Annotate standardized file
sbatch scripts/annotateVCF.sh original_vcfs/pe_pop3_filt.standardized.vcf.gz pe_pop3_filt
```

### 3. Calculate VCF Statistics
Generate quality metrics using bcftools and vcftools. The arguments it takes is path/to/vcf.

```bash
# Statistics for files being compared
sbatch scripts/getSNPstats.sh vcfs/pe_pop3.annotated.vcf.gz
sbatch scripts/getSNPstats.sh vcfs/populations.all.annotated.vcf.gz
sbatch scripts/getSNPstats.sh vcfs/pe_pop3_filt.annotated.vcf.gz
```

### 4. Visualize QC Statistics
Generate plots to compare filtering steps (skip invariant sites vcf).

```bash
# Plot stats for all files except populations.all
for vcf in vcfs/*.gz; do 
  if [[ "$vcf" != "vcfs/populations.all.annotated.vcf.gz" ]]; then
    sbatch scripts/plotRawStats.sh $vcf
  fi
done
```

### 5. Final Filtering
Identifies and removes samples with >30% missingness, re-annotates the cleaned VCF, applies individual genotype filters (depth <3, quality <20), and filters sites with >20% missingness or MAF ≤0.03. The arguments it takes is path/to/vcf.

```bash
sbatch scripts/02-final-filter.sh vcfs/new.filt.vcf.gz
```

### 6. Population Structure Analysis
Generate PCA and admixture plots. The arguments it takes is path/to/vcf, then window size, step size and then linkage correlation (r2-threshold).

```bash
# Arguments: <vcf> <max_PC> <max_K> <maf_threshold>
sbatch scripts/getPopStruct.sh vcfs/new.filt.final.vcf.gz 50 10 0.2
```

### 7. Merge with Invariant Sites
Create final VCF with both variant and invariant sites to be pixy/piawka compliant for calculations.

```bash
sbatch scripts/mergeInvariants.sh
```

### 8. Population Genetics Statistics
Maps samples to populations, then calculates pairwise FST, nucleotide diversity (π), and Dxy statistics using piawka with a maximum of 30% missing data per site per population group. The arguments it takes is path/to/vcf.

```bash
sbatch scripts/runPiawka.sh vcfs/final.merged.vcf.gz
```

### 9. Microbiome Analysis (TBA)
Process amplicon data and correlate with host genetics.

## Main Tools
mamba environments in .yml files can be found in `mamba_environments/`. 
| Name | Version | Channel |
|------|---------|---------|
| ADMIXTURE | 1.3.0 | dardel_module |
| bcftools | 1.20 | dardel_module |
| piawka | 0.8.11 | bioconda |
| plink | 2.00a5.14 | dardel_module |
| r-base | 4.3-4.5.1 | conda-forge |
| vcftools | 0.1.16 | dardel_module |