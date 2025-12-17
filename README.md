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
#### VCFs of Importance
```
VCF                             	 samples	 records 	 invariants	 SNPs 	 indels	 others
# For population genomics	
final.merged.vcf.gz             	      93	 18056199	   18035778	 20421	      0	      0	

# Newly created with updated filtering
new.filt.final.vcf.gz           	      93	    20421	          0	 20421	      0	      0	

# Past final version
pe_pop3_filt.annotated.vcf.gz   	      96	      174	          0	   174	      0	      0	

# Raw VCF with only biallelic snps 
pe_pop3.annotated.vcf.gz        	      96	    98825	          0	 98825	      0	      0	

# Raw VCF with biallelic snps + invariants
populations.all.annotated.vcf.gz	      96	 18134603	   18035778	 98825	      0	      0	

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

## Working with the Data Frames
In `02-master-data-frame.R` there are two R objects of importance: `master_frames` and `df_all`
_example_

```R
# All mantel dataframes stored in list at their taxonomic level.
master_frames <- master_list

# Save the object to an .Rds file 
saveRDS(master_frames, file = "data/master_frames.Rds") 

# Combine all taxonomic levels into one dataframe
df_all <- bind_rows(
  master_frames$Phylum %>% mutate(taxon_level = "Phylum"),
  master_frames$Class %>% mutate(taxon_level = "Class"),
  master_frames$Order %>% mutate(taxon_level = "Order"),
  master_frames$Family %>% mutate(taxon_level = "Family"),
  master_frames$Genus %>% mutate(taxon_level = "Genus"),
  master_frames$Species %>% mutate(taxon_level = "Species")
)

# Save as .csv and .Rds
write_csv(df_all, "data/mantel_all_taxa.csv")
saveRDS(df_all, file = "data/mantel_all_taxa.Rds")
```

These can be loaded using this code
```R
# You can load it back with: 
master_frames <- readRDS("data/master_frames.Rds")
df_all <- readRDS("data/mantel_all_taxa.Rds")
```

Preview of `mantel_all_taxa.Rds`
```R
> # Show dimensions
> dim(df_all)
[1] 384   9

> # Show number of samples per taxa
> df_all %>% count(taxon_level)
# A tibble: 6 × 2
  taxon_level     n
  <chr>       <int>
1 Class          64
2 Family         64
3 Genus          64
4 Order          64
5 Phylum         64
6 Species        64

> # Shows three rows from each taxonomic level
> df_all %>% group_by(taxon_level) %>% slice_head(n = 3)
# A tibble: 18 × 9
# Groups:   taxon_level [6]
   pair  pop1.id pop2.id   fst mean_bray.AMF mean_bray.ITS mean_bray.16S km.diff taxon_level
   <chr> <chr>   <chr>   <dbl>         <dbl>         <dbl>         <dbl>   <dbl> <chr>      
 1 B1-B1 B1      B1      0            0.0389        0.328          0.240     0   Class      
 2 B1-B2 B1      B2      0.258        0.0505        0.363          0.268    31.2 Class      
 3 B1-F1 B1      F1      0.544        0.0366        0.357          0.332     1.9 Class      
 4 B1-B1 B1      B1      0            0.355         0.540          0.445     0   Family     
 5 B1-B2 B1      B2      0.258        0.416         0.646          0.486    31.2 Family     
 6 B1-F1 B1      F1      0.544        0.468         0.721          0.593     1.9 Family     
 7 B1-B1 B1      B1      0            0.356         0.596          0.548     0   Genus      
 8 B1-B2 B1      B2      0.258        0.416         0.745          0.600    31.2 Genus      
 9 B1-F1 B1      F1      0.544        0.468         0.774          0.719     1.9 Genus      
10 B1-B1 B1      B1      0            0.237         0.488          0.375     0   Order      
11 B1-B2 B1      B2      0.258        0.344         0.545          0.416    31.2 Order      
12 B1-F1 B1      F1      0.544        0.403         0.543          0.493     1.9 Order      
13 B1-B1 B1      B1      0            0             0.0637         0.191     0   Phylum     
14 B1-B2 B1      B2      0.258        0             0.0784         0.230    31.2 Phylum     
15 B1-F1 B1      F1      0.544        0             0.114          0.250     1.9 Phylum     
16 B1-B1 B1      B1      0            0.370         0.640          0.570     0   Species    
17 B1-B2 B1      B2      0.258        0.427         0.806          0.632    31.2 Species    
18 B1-F1 B1      F1      0.544        0.477         0.872          0.748     1.9 Species    
```

Preview of `master_frames.Rds` and how to use it.
```R
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

...

$Species
# A tibble: 64 × 8
   pair  pop1.id pop2.id   fst mean_bray.AMF mean_bray.ITS mean_bray.16S km.diff
   <chr> <chr>   <chr>   <dbl>         <dbl>         <dbl>         <dbl>   <dbl>
 1 B1-B1 B1      B1      0             0.370         0.640         0.570     0  
 2 B1-B2 B1      B2      0.258         0.427         0.806         0.632    31.2
 3 B1-F1 B1      F1      0.544         0.477         0.872         0.748     1.9
 4 B1-F2 B1      F2      0.261         0.491         0.836         0.675     2.3
 5 B1-R1 B1      R1      0.602         0.457         0.837         0.806     9.2
 6 B1-R2 B1      R2      0.479         0.542         0.819         0.685    24.2
 7 B1-U1 B1      U1      0.402         0.458         0.837         0.662    10.2
 8 B1-U2 B1      U2      0.409         0.409         0.830         0.647    11.9
 9 B2-B1 B2      B1      0.258         0.427         0.806         0.632    31.2
10 B2-B2 B2      B2      0             0.225         0.669         0.530     0  
# ℹ 54 more rows
# ℹ Use `print(n = ...)` to see more rows

> # Subset to specific taxonomic level: Order
> df_order <- master_frames$Order

> dim(df_order)
[1] 64  8

> df_order
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

## Repo size
Note: `processDADA2` contains sequence data and is compressed. Original data ~50GB
```
du -sh .
14G     .

du -sh *
204K    01-pruneLD
1.6M    02-PCA
4.9M    03-ADMIXTURE
16K     04-piawka
134M    05-mantel-overview
24K     README.md
6.6M    data
1.4M    doc
8.0K    fst_dxy_matrix.xlsx
20K     mamba_environments
250M    original_vcfs
56K     plots
11G     processDADA2
100K    scripts
14M     slurm-logs
1.7G    stats
686M    vcfs
```

## File Structure 
```
.
├── 01-pruneLD
│   └── new.filt.final-50-10-0.2
│       ├── new.filt.final-50-10-0.2.log
│       ├── new.filt.final-50-10-0.2.prune.in
│       └── new.filt.final-50-10-0.2.prune.out
├── 02-PCA
│   └── new.filt.final-50-10-0.2
│       ├── PCA_1vs2.png
│       ├── PCA_3vs4.png
│       ├── PCA_5vs6.png
│       ├── PCA_scree.png
│       ├── new.filt.final-50-10-0.2.bed
│       ├── new.filt.final-50-10-0.2.bim
│       ├── new.filt.final-50-10-0.2.eigenval
│       ├── new.filt.final-50-10-0.2.eigenvec
│       ├── new.filt.final-50-10-0.2.fam
│       ├── new.filt.final-50-10-0.2.human-chr.bim
│       ├── new.filt.final-50-10-0.2.log
│       └── samples.txt
├── 03-ADMIXTURE
│   └── new.filt.final-50-10-0.2
│       ├── admixture-2.png
│       ├── admixture-3.png
│       ├── admixture-4.png
│       ├── admixture-5.png
│       ├── admixture-6.png
│       ├── admixture-7.png
│       ├── admixture-8.png
│       ├── admixture-9.png
│       ├── new.filt.final-50-10-0.2-log1.out
│       ├── new.filt.final-50-10-0.2-log2.out
│       ├── new.filt.final-50-10-0.2-log3.out
│       ├── new.filt.final-50-10-0.2-log4.out
│       ├── new.filt.final-50-10-0.2-log5.out
│       ├── new.filt.final-50-10-0.2-log6.out
│       ├── new.filt.final-50-10-0.2-log7.out
│       ├── new.filt.final-50-10-0.2-log8.out
│       ├── new.filt.final-50-10-0.2-log9.out
│       ├── new.filt.final-50-10-0.2.1.P
│       ├── new.filt.final-50-10-0.2.1.Q
│       ├── new.filt.final-50-10-0.2.2.P
│       ├── new.filt.final-50-10-0.2.2.Q
│       ├── new.filt.final-50-10-0.2.3.P
│       ├── new.filt.final-50-10-0.2.3.Q
│       ├── new.filt.final-50-10-0.2.4.P
│       ├── new.filt.final-50-10-0.2.4.Q
│       ├── new.filt.final-50-10-0.2.5.P
│       ├── new.filt.final-50-10-0.2.5.Q
│       ├── new.filt.final-50-10-0.2.6.P
│       ├── new.filt.final-50-10-0.2.6.Q
│       ├── new.filt.final-50-10-0.2.7.P
│       ├── new.filt.final-50-10-0.2.7.Q
│       ├── new.filt.final-50-10-0.2.8.P
│       ├── new.filt.final-50-10-0.2.8.Q
│       ├── new.filt.final-50-10-0.2.9.P
│       └── new.filt.final-50-10-0.2.9.Q
├── 04-piawka
│   ├── final.merged-CnN1qz.groups
│   ├── final.merged-tRowBE.samples
│   └── final.merged.tsv
├── 05-mantel-overview
│   ├── mantel-16S
│   │   ├── 01-16S-fst-bray-no-self-Phylum.jpg
│   │   ├── 01-16S-fst-bray-no-self-Phylum.svg
│   │   ├── 02-16S-fst-bray-no-self-Class.jpg
│   │   ├── 02-16S-fst-bray-no-self-Class.svg
│   │   ├── 03-16S-fst-bray-no-self-Order.jpg
│   │   ├── 03-16S-fst-bray-no-self-Order.svg
│   │   ├── 04-16S-fst-bray-no-self-Family.jpg
│   │   ├── 04-16S-fst-bray-no-self-Family.svg
│   │   ├── 05-16S-fst-bray-no-self-Genus.jpg
│   │   ├── 05-16S-fst-bray-no-self-Genus.svg
│   │   ├── 06-16S-fst-bray-no-self-Species.jpg
│   │   └── 06-16S-fst-bray-no-self-Species.svg
│   ├── mantel-AMF
│   │   ├── 01-AMF-fst-bray-no-self-Phylum.jpg
│   │   ├── 01-AMF-fst-bray-no-self-Phylum.svg
│   │   ├── 02-AMF-fst-bray-no-self-Class.jpg
│   │   ├── 02-AMF-fst-bray-no-self-Class.svg
│   │   ├── 03-AMF-fst-bray-no-self-Order.jpg
│   │   ├── 03-AMF-fst-bray-no-self-Order.svg
│   │   ├── 04-AMF-fst-bray-no-self-Family.jpg
│   │   ├── 04-AMF-fst-bray-no-self-Family.svg
│   │   ├── 05-AMF-fst-bray-no-self-Genus.jpg
│   │   ├── 05-AMF-fst-bray-no-self-Genus.svg
│   │   ├── 06-AMF-fst-bray-no-self-Species.jpg
│   │   └── 06-AMF-fst-bray-no-self-Species.svg
│   ├── mantel-ITS
│   │   ├── 01-ITS-fst-bray-no-self-Phylum.jpg
│   │   ├── 01-ITS-fst-bray-no-self-Phylum.svg
│   │   ├── 02-ITS-fst-bray-no-self-Class.jpg
│   │   ├── 02-ITS-fst-bray-no-self-Class.svg
│   │   ├── 03-ITS-fst-bray-no-self-Order.jpg
│   │   ├── 03-ITS-fst-bray-no-self-Order.svg
│   │   ├── 04-ITS-fst-bray-no-self-Family.jpg
│   │   ├── 04-ITS-fst-bray-no-self-Family.svg
│   │   ├── 05-ITS-fst-bray-no-self-Genus.jpg
│   │   ├── 05-ITS-fst-bray-no-self-Genus.svg
│   │   ├── 06-ITS-fst-bray-no-self-Species.jpg
│   │   └── 06-ITS-fst-bray-no-self-Species.svg
│   ├── mantel-together
│   │   ├── 01-mantel-comparison-Phylum.jpg
│   │   ├── 01-mantel-comparison-Phylum.svg
│   │   ├── 02-mantel-comparison-Class.jpg
│   │   ├── 02-mantel-comparison-Class.svg
│   │   ├── 03-mantel-comparison-Order.jpg
│   │   ├── 03-mantel-comparison-Order.svg
│   │   ├── 04-mantel-comparison-Family.jpg
│   │   ├── 04-mantel-comparison-Family.svg
│   │   ├── 05-mantel-comparison-Genus.jpg
│   │   ├── 05-mantel-comparison-Genus.svg
│   │   ├── 06-mantel-comparison-Species.jpg
│   │   └── 06-mantel-comparison-Species.svg
│   ├── mantel_summary.xlsx
│   ├── ridge-combo-16S
│   │   ├── 01-16S-fst-bray-ridge-Phylum.jpg
│   │   ├── 01-16S-fst-bray-ridge-Phylum.svg
│   │   ├── 02-16S-fst-bray-ridge-Class.jpg
│   │   ├── 02-16S-fst-bray-ridge-Class.svg
│   │   ├── 03-16S-fst-bray-ridge-Order.jpg
│   │   ├── 03-16S-fst-bray-ridge-Order.svg
│   │   ├── 04-16S-fst-bray-ridge-Family.jpg
│   │   ├── 04-16S-fst-bray-ridge-Family.svg
│   │   ├── 05-16S-fst-bray-ridge-Genus.jpg
│   │   ├── 05-16S-fst-bray-ridge-Genus.svg
│   │   ├── 06-16S-fst-bray-ridge-Species.jpg
│   │   └── 06-16S-fst-bray-ridge-Species.svg
│   ├── ridge-combo-AMF
│   │   ├── 01-AMF-fst-bray-ridge-Phylum.jpg
│   │   ├── 01-AMF-fst-bray-ridge-Phylum.svg
│   │   ├── 02-AMF-fst-bray-ridge-Class.jpg
│   │   ├── 02-AMF-fst-bray-ridge-Class.svg
│   │   ├── 03-AMF-fst-bray-ridge-Order.jpg
│   │   ├── 03-AMF-fst-bray-ridge-Order.svg
│   │   ├── 04-AMF-fst-bray-ridge-Family.jpg
│   │   ├── 04-AMF-fst-bray-ridge-Family.svg
│   │   ├── 05-AMF-fst-bray-ridge-Genus.jpg
│   │   ├── 05-AMF-fst-bray-ridge-Genus.svg
│   │   ├── 06-AMF-fst-bray-ridge-Species.jpg
│   │   └── 06-AMF-fst-bray-ridge-Species.svg
│   └── ridge-combo-ITS
│       ├── 01-ITS-fst-bray-ridge-Phylum.jpg
│       ├── 01-ITS-fst-bray-ridge-Phylum.svg
│       ├── 02-ITS-fst-bray-ridge-Class.jpg
│       ├── 02-ITS-fst-bray-ridge-Class.svg
│       ├── 03-ITS-fst-bray-ridge-Order.jpg
│       ├── 03-ITS-fst-bray-ridge-Order.svg
│       ├── 04-ITS-fst-bray-ridge-Family.jpg
│       ├── 04-ITS-fst-bray-ridge-Family.svg
│       ├── 05-ITS-fst-bray-ridge-Genus.jpg
│       ├── 05-ITS-fst-bray-ridge-Genus.svg
│       ├── 06-ITS-fst-bray-ridge-Species.jpg
│       └── 06-ITS-fst-bray-ridge-Species.svg
├── README.md
├── data
│   ├── 16S_abundance
│   │   ├── Class_abund.csv
│   │   ├── Family_abund.csv
│   │   ├── Genus_abund.csv
│   │   ├── Order_abund.csv
│   │   ├── Phylum_abund.csv
│   │   └── Species_abund.csv
│   ├── AMF_abundance
│   │   ├── Class_abund.tsv
│   │   ├── Family_abund.tsv
│   │   ├── Genus_abund.tsv
│   │   ├── Order_abund.tsv
│   │   ├── Phylum_abund.tsv
│   │   └── Species_abund.tsv
│   ├── ITS_abundance
│   │   ├── Class_abund.csv
│   │   ├── Family_abund.csv
│   │   ├── Genus_abund.csv
│   │   ├── Order_abund.csv
│   │   ├── Phylum_abund.csv
│   │   └── Species_abund.csv
│   ├── mantel_class.csv
│   ├── mantel_family.csv
│   ├── mantel_genus.csv
│   ├── mantel_order.csv
│   ├── mantel_phylum.csv
│   └── mantel_species.csv
├── doc
│   ├── S_dulcamara_popmap_1.txt
│   ├── S_dulcamara_sample_info_1.txt
│   ├── chr_map_replace-invariants.txt
│   ├── geo_distance.tsv
│   ├── meta-microeco.tsv
│   ├── metadata-AMF.tsv
│   ├── pop_map.tsv
│   └── uniq_pop_meta.tsv
├── fst_dxy_matrix.xlsx
├── mamba_environments
│   ├── piawka-environment.yml
│   ├── r-arena-environment.yml
│   └── reproduce_dada2-environment.yml
├── original_vcfs
│   ├── pe_pop3.snps.vcf.gz
│   ├── pe_pop3.snps.vcf.gz.tbi
│   ├── pe_pop3_filt.snps.vcf.gz
│   ├── pe_pop3_filt.standardized.vcf.gz
│   ├── pe_pop3_filt.standardized.vcf.gz.tbi
│   ├── populations.all.vcf.gz
│   └── populations.all.vcf.gz.tbi
├── plots
│   ├── dxy_small_heatmap.png
│   └── fst_small_heatmap.png
├── processDADA2
│   ├── amplicon_data.tar.gz
│   ├── amplicon_rds.tar.gz
│   ├── scripts
│   │   ├── r-scripts
│   │   │   ├── create16Sotu.R
│   │   │   ├── createITSotu.R
│   │   │   ├── procesSDADA2-16S.R
│   │   │   └── processDADA2-ITS.R
│   │   ├── runDADA2-16S.sh
│   │   └── runDADA2-ITS.sh
│   └── tax
│       ├── sh_general_release_dynamic_s_all_10.05.2021.fasta
│       ├── silva_nr99_v138.1_train_set.fa.gz
│       └── silva_species_assignment_v138.1.fa.gz
├── scripts
│   ├── 01-initial-filter.sh
│   ├── 02-final-filter.sh
│   ├── annotateVCF.sh
│   ├── getBCFstats.sh
│   ├── getPopStruct.sh
│   ├── getSNPstats.sh
│   ├── mergeInvariants.sh
│   ├── plotRawStats.sh
│   ├── r-scripts
│   │   ├── 01_mantel.R
│   │   ├── 02-master-data-frame.R
│   │   ├── fst_dxy_hc.R
│   │   ├── master-data-frame-functions.R
│   │   ├── plotADMIXTURE.R
│   │   ├── plotPCA.R
│   │   ├── plotRawStats.R
│   │   └── plotVCFTOOLSstats.R
│   ├── runPiawka.sh
│   └── tarball-it.sh
├── slurm-logs
│   ├── PopStruct
│   ├── R
│   ├── annotate
│   ├── piawka
│   ├── stats
│   ├── tarball
│   └── wrangle
├── stats
│   ├── new.filt
│   │   ├── new.filt.frq
│   │   ├── new.filt.idepth
│   │   ├── new.filt.imiss
│   │   ├── new.filt.lmiss
│   │   ├── plots
│   │   │   ├── AC_density.png
│   │   │   ├── AF_density.png
│   │   │   ├── AN_density.png
│   │   │   ├── ExcHet_density.png
│   │   │   ├── F_MISSING.png
│   │   │   ├── idepth_boxplot.png
│   │   │   ├── imiss_boxplot.png
│   │   │   ├── lmiss_density.png
│   │   │   ├── maf_density.png
│   │   │   └── maf_threshold_scree.png
│   │   └── raw_stats-new.filt.tsv
│   ├── new.filt.clean.anno
│   │   ├── new.filt.clean.anno.frq
│   │   ├── new.filt.clean.anno.idepth
│   │   ├── new.filt.clean.anno.imiss
│   │   ├── new.filt.clean.anno.lmiss
│   │   └── raw_stats-new.filt.clean.anno.tsv
│   ├── new.filt.final
│   │   ├── new.filt.final.frq
│   │   ├── new.filt.final.idepth
│   │   ├── new.filt.final.imiss
│   │   ├── new.filt.final.lmiss
│   │   └── raw_stats-new.filt.final.tsv
│   ├── new.raw.anno
│   │   ├── new.raw.anno.frq
│   │   ├── new.raw.anno.idepth
│   │   ├── new.raw.anno.imiss
│   │   ├── new.raw.anno.lmiss
│   │   ├── plots
│   │   │   ├── AC_density.png
│   │   │   ├── AF_density.png
│   │   │   ├── AN_density.png
│   │   │   ├── ExcHet_density.png
│   │   │   ├── F_MISSING.png
│   │   │   ├── idepth_boxplot.png
│   │   │   ├── imiss_boxplot.png
│   │   │   ├── lmiss_density.png
│   │   │   ├── maf_density.png
│   │   │   └── maf_threshold_scree.png
│   │   └── raw_stats-new.raw.anno.tsv
│   ├── pe_pop3.annotated
│   │   ├── pe_pop3.annotated.frq
│   │   ├── pe_pop3.annotated.idepth
│   │   ├── pe_pop3.annotated.imiss
│   │   ├── pe_pop3.annotated.lmiss
│   │   ├── plots
│   │   │   ├── AC_density.png
│   │   │   ├── AF_density.png
│   │   │   ├── AN_density.png
│   │   │   ├── ExcHet_density.png
│   │   │   ├── F_MISSING.png
│   │   │   ├── idepth_boxplot.png
│   │   │   ├── imiss_boxplot.png
│   │   │   ├── lmiss_density.png
│   │   │   ├── maf_density.png
│   │   │   └── maf_threshold_scree.png
│   │   └── raw_stats-pe_pop3.annotated.tsv
│   ├── pe_pop3_filt.annotated
│   │   ├── pe_pop3_filt.annotated.frq
│   │   ├── pe_pop3_filt.annotated.idepth
│   │   ├── pe_pop3_filt.annotated.imiss
│   │   ├── pe_pop3_filt.annotated.lmiss
│   │   ├── plots
│   │   │   ├── AC_density.png
│   │   │   ├── AF_density.png
│   │   │   ├── AN_density.png
│   │   │   ├── ExcHet_density.png
│   │   │   ├── F_MISSING.png
│   │   │   ├── idepth_boxplot.png
│   │   │   ├── imiss_boxplot.png
│   │   │   ├── lmiss_density.png
│   │   │   ├── maf_density.png
│   │   │   └── maf_threshold_scree.png
│   │   └── raw_stats-pe_pop3_filt.annotated.tsv
│   ├── populations.all.annotated
│   │   ├── plots
│   │   │   ├── AC_density.png
│   │   │   └── AN_density.png
│   │   ├── populations.all.annotated.frq
│   │   ├── populations.all.annotated.idepth
│   │   ├── populations.all.annotated.imiss
│   │   ├── populations.all.annotated.lmiss
│   │   └── raw_stats-populations.all.annotated.tsv
│   ├── samples_to_remove.txt
│   └── vcf_stats.tsv
└── vcfs
    ├── final.invariants.vcf.gz
    ├── final.invariants.vcf.gz.tbi
    ├── final.merged.vcf.gz
    ├── final.merged.vcf.gz.tbi
    ├── new.filt.clean.anno.vcf.gz
    ├── new.filt.clean.anno.vcf.gz.tbi
    ├── new.filt.final.vcf.gz
    ├── new.filt.final.vcf.gz.tbi
    ├── new.filt.vcf.gz
    ├── new.filt.vcf.gz.tbi
    ├── new.raw.anno.vcf.gz
    ├── new.raw.anno.vcf.gz.tbi
    ├── pe_pop3.annotated.vcf.gz
    ├── pe_pop3.annotated.vcf.gz.tbi
    ├── pe_pop3_filt.annotated.vcf.gz
    ├── pe_pop3_filt.annotated.vcf.gz.tbi
    ├── populations.all.annotated.vcf.gz
    └── populations.all.annotated.vcf.gz.tbi
```
