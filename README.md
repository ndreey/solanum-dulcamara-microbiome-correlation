# Solanum dulcamara Population Genomics and Microbiome Analysis

This project investigates the relationship between host plant population genetics and root-associated microbiome composition in *Solanum dulcamara* across natural populations. Using RADseq data and 16S/ITS amplicon sequencing, we test whether host genetic distance (FST) correlates with microbiome beta-diversity (Bray-Curtis dissimilarity) through Mantel correlation tests. The analysis pipeline includes vcf filtering, population structure analysis (PCA, admixture), calculation of population differentiation (FST), nucleotide divergence (Dxy), and microbiome community. We compare bacterial (16S), fungal (ITS), and arbuscular mycorrhizal fungi (AMF) communities across taxonomic levels to determine which microbial groups show the strongest correlation with host genetic structure.

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

### 9. Microbiome Analysis
Process amplicon data and correlate with host genetics and distances.
Used R packages:
```
tidyverse     2.0.0
vegan         2.8-0
scales        1.4.0
ggpubr        0.6.2 
ggridges      0.5.7
patchwork     1.3.2
writexl       1.5.4
pheatmap      1.0.13
paletteer     1.6.0
```
Run these scripts within Rstudio:
- `01-mantel.R`: This will create multiple plots such as ridge plots, scatter plots to better visualize the mean Bray-Curtis data.
- `02-master-data-frame.R`: This will create a list where each element is a dataframe for each taxonomic level. Holding FST, km.diff, and mean Bray-Curtis distances between each pairwise comparison.

_example_
```
> master_frames
$Phylum
# A tibble: 64 × 8
   pair  pop1.id pop2.id   fst mean_bray.AMF mean_bray.ITS mean_bray.16S km.diff
   <chr> <chr>   <chr>   <dbl>         <dbl>         <dbl>         <dbl>   <dbl>
 1 B1-B1 B1      B1      0                 0        0.0637         0.191     0  
 2 B1-B2 B1      B2      0.258             0        0.0784         0.230    31.2
 3 B1-F1 B1      F1      0.544             0        0.114          0.250     1.9
 4 B1-F2 B1      F2      0.261             0        0.0811         0.201     2.3
 5 B1-R1 B1      R1      0.602             0        0.0850         0.345     9.2
 6 B1-R2 B1      R2      0.479             0        0.0759         0.201    24.2
 7 B1-U1 B1      U1      0.402             0        0.0615         0.183    10.2
 8 B1-U2 B1      U2      0.409             0        0.0765         0.211    11.9
 9 B2-B1 B2      B1      0.258             0        0.0784         0.230    31.2
10 B2-B2 B2      B2      0                 0        0.0810         0.200     0  
# ℹ 54 more rows
# ℹ Use `print(n = ...)` to see more rows

$Class
# A tibble: 64 × 8
   pair  pop1.id pop2.id   fst mean_bray.AMF mean_bray.ITS mean_bray.16S km.diff
   <chr> <chr>   <chr>   <dbl>         <dbl>         <dbl>         <dbl>   <dbl>
 1 B1-B1 B1      B1      0            0.0389         0.328         0.240     0  
 2 B1-B2 B1      B2      0.258        0.0505         0.363         0.268    31.2
 3 B1-F1 B1      F1      0.544        0.0366         0.357         0.332     1.9
 4 B1-F2 B1      F2      0.261        0.0377         0.305         0.272     2.3
 5 B1-R1 B1      R1      0.602        0.0493         0.318         0.437     9.2
 6 B1-R2 B1      R2      0.479        0.0641         0.291         0.268    24.2
 7 B1-U1 B1      U1      0.402        0.0533         0.371         0.251    10.2
 8 B1-U2 B1      U2      0.409        0.0446         0.345         0.264    11.9
 9 B2-B1 B2      B1      0.258        0.0505         0.363         0.268    31.2
10 B2-B2 B2      B2      0            0.0490         0.304         0.235     0  
# ℹ 54 more rows
# ℹ Use `print(n = ...)` to see more rows

$Order
# A tibble: 64 × 8
   pair  pop1.id pop2.id   fst mean_bray.AMF mean_bray.ITS mean_bray.16S km.diff
   <chr> <chr>   <chr>   <dbl>         <dbl>         <dbl>         <dbl>   <dbl>
 1 B1-B1 B1      B1      0             0.237         0.488         0.375     0  
 2 B1-B2 B1      B2      0.258         0.344         0.545         0.416    31.2
 3 B1-F1 B1      F1      0.544         0.403         0.543         0.493     1.9
 4 B1-F2 B1      F2      0.261         0.411         0.495         0.450     2.3
 5 B1-R1 B1      R1      0.602         0.387         0.515         0.587     9.2
 6 B1-R2 B1      R2      0.479         0.502         0.456         0.438    24.2
 7 B1-U1 B1      U1      0.402         0.356         0.570         0.419    10.2
 8 B1-U2 B1      U2      0.409         0.301         0.481         0.413    11.9
 9 B2-B1 B2      B1      0.258         0.344         0.545         0.416    31.2
10 B2-B2 B2      B2      0             0.202         0.411         0.350     0  
# ℹ 54 more rows
# ℹ Use `print(n = ...)` to see more rows
```

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

## Reproduce DADA2 tables.
For this, one will have to move the cwd to the `processDADA2/` directory
and run scripts from there.