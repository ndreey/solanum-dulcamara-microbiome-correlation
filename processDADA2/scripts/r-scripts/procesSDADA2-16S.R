library(dada2)
library(RcppParallel)

# Set threads
setThreadOptions(numThreads = 8)

# Create output directory if it doesn't exist
dir.create("amplicon_rds/16S/", showWarnings = FALSE, recursive = TRUE)

# Path to demultiplexed reads
path_raw <- "amplicon_data/16S/raw"
path_filt <- "amplicon_data/16S/filt"

message("Getting files")
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fq and SAMPLENAME_2.fq
fnFs <- sort(list.files(path_raw, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_2.fq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

message("Setting names")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_filt, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path_filt, paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# ============ FILTER AND TRIM ============
if (file.exists("amplicon_rds/16S/out.rds")) {
  message("Loading existing filter results from RDS")
  out <- readRDS("amplicon_rds/16S/out.rds")
} else {
  message("Standard filtering starting")
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,210),
                maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                compress=TRUE, multithread=8)
  saveRDS(out, "amplicon_rds/16S/out.rds")
  message("Filter complete and saved")
}
print(dim(out))

# ============ ERROR MODELS ============
if (file.exists("amplicon_rds/16S/errF.rds")) {
  message("Loading existing forward error model from RDS")
  errF <- readRDS("amplicon_rds/16S/errF.rds")
} else {
  message("Learning forward error model")
  errF <- learnErrors(filtFs, multithread=8)
  saveRDS(errF, "amplicon_rds/16S/errF.rds")
  message("Forward error model saved")
}

if (file.exists("amplicon_rds/16S/errR.rds")) {
  message("Loading existing reverse error model from RDS")
  errR <- readRDS("amplicon_rds/16S/errR.rds")
} else {
  message("Learning reverse error model")
  errR <- learnErrors(filtRs, multithread=8)
  saveRDS(errR, "amplicon_rds/16S/errR.rds")
  message("Reverse error model saved")
}

# ============ SAMPLE INFERENCE ============
if (file.exists("amplicon_rds/16S/dadaFs.rds") && file.exists("amplicon_rds/16S/dadaRs.rds")) {
  message("Loading existing dada results from RDS")
  dadaFs <- readRDS("amplicon_rds/16S/dadaFs.rds")
  dadaRs <- readRDS("amplicon_rds/16S/dadaRs.rds")
} else {
  message("Applying core sample inference")
  dadaFs <- dada(filtFs, err=errF, multithread=8)
  saveRDS(dadaFs, "amplicon_rds/16S/dadaFs.rds")
  message("Sample inference complete and saved for forward")
  dadaRs <- dada(filtRs, err=errR, multithread=8)
  saveRDS(dadaRs, "amplicon_rds/16S/dadaRs.rds")
  message("Sample inference complete and saved for reverse")

}

# ============ MERGE READS ============
if (file.exists("amplicon_rds/16S/mergers.rds")) {
  message("Loading existing merged reads from RDS")
  mergers <- readRDS("amplicon_rds/16S/mergers.rds")
} else {
  message("Merging reads")
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  saveRDS(mergers, "amplicon_rds/16S/mergers.rds")
  message("Merge complete and saved")
}

# ============ SEQUENCE TABLE ============
if (file.exists("amplicon_rds/16S/seqtab.rds")) {
  message("Loading existing sequence table from RDS")
  seqtab <- readRDS("amplicon_rds/16S/seqtab.rds")
} else {
  message("Creating sequence table")
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, "amplicon_rds/16S/seqtab.rds")
  message("Sequence table saved")
}
message(paste("Dimensions:", dim(seqtab)[1], "samples,", dim(seqtab)[2], "ASVs"))
message("Sequence length distribution:")
print(table(nchar(getSequences(seqtab))))

# ============ REMOVE CHIMERAS ============
if (file.exists("amplicon_rds/16S/seqtab.nochim.rds")) {
  message("Loading existing chimera-free sequence table from RDS")
  seqtab.nochim <- readRDS("amplicon_rds/16S/seqtab.nochim.rds")
} else {
  message("Removing chimeras")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
  saveRDS(seqtab.nochim, "amplicon_rds/16S/seqtab.nochim.rds")
  message("Chimera removal complete and saved")
}
message(paste("Dimensions after chimera removal:", dim(seqtab.nochim)[1], "samples,", dim(seqtab.nochim)[2], "ASVs"))
message(paste("Proportion of reads retained:", round(sum(seqtab.nochim)/sum(seqtab), 3)))

# ============ TRACK READS ============
if (file.exists("amplicon_rds/16S/track.rds")) {
  message("Loading existing read tracking from RDS")
  track <- readRDS("amplicon_rds/16S/track.rds")
} else {
  message("Track reads through pipeline")
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
                 sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  saveRDS(track, "amplicon_rds/16S/track.rds")
  message("Read tracking saved")
}
print(head(track))

# ============ ASSIGN TAXONOMY ============
if (file.exists("amplicon_rds/16S/taxa.rds")) {
  message("Loading existing taxonomy from RDS")
  taxa <- readRDS("amplicon_rds/16S/taxa.rds")
} else {
  message("Assign taxonomy")
  taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread=8)
  saveRDS(taxa, "amplicon_rds/16S/taxa.rds")
  message("Taxonomy assignment complete and saved")
}

# Check if species assignment already done
if (file.exists("amplicon_rds/16S/taxa_species.rds")) {
  message("Loading existing species-level taxonomy from RDS")
  taxa_species <- readRDS("amplicon_rds/16S/taxa_species.rds")
} else {
  message("Assign species")
  taxa_species <- addSpecies(taxa, "tax/silva_species_assignment_v138.1.fa.gz")
  saveRDS(taxa_species, "amplicon_rds/16S/taxa_species.rds")
  message("Species assignment complete and saved")
}

message("Checking taxa")
taxa.print <- taxa_species  # Use the version with species
rownames(taxa.print) <- NULL
print(head(taxa.print))

message("Checking taxa")
taxa.print <- taxa
rownames(taxa.print) <- NULL
print(head(taxa.print))

message("DADA2 pipeline complete! All objects saved to amplicon_rds/16S/")