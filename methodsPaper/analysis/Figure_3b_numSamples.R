## Script to be run on HPC with access to per-sample peak sets and group peak sets

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

## Functions ## 
ReadPeakstoGRanges <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue")
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end), signalValue = df$signalValue,pValue = df$pValue,qValue = df$qValue )
  return(gr)
}

## Directories
subsetDir <- "/path/to/output/file"
samplesPeaksDir <- "/path/to/samples/peaks"
metaDir <- "/path/to/list/samples/"
cf <- "selected cell-fraction"
groupPeaksDir <- "/path/to/group/peaks"

# Gather samples peak sets
peak_files <- list.files(path = samplesPeaksDir, pattern = "(_merged)?\\.narrowPeak\\.filt$", full.names = TRUE)
length(peak_files)

# Select samples for group analysis
samples <- read.table(file.path(paste0(metaDir,"downAnalysisSamples_",cf,".txt")))[,1]
length(samples)

# Extract file names without extensions for labeling
file_labels <- basename(peak_files)
file_labels <- gsub(".narrowPeak(\\.filt)?$", "", file_labels)
names(peak_files) <- file_labels
# Check peak files for all samples in list exist
peak_files <- peak_files[names(peak_files) %in% samples]
print("Names of samples peaks found: ")
length(peak_files)
setdiff(samples,names(peak_files))

# Read peak files as genomic ranges
peak_sets <- lapply(peak_files, ReadPeakstoGRanges)

# Read cell-group peaks
peak_file_group <-paste0(groupPeaksDir,cf,"_5e-2.pvalsorted.narrowPeak.filt")
peak_set_group <- ReadPeakstoGRanges(peak_file_group)
names(peak_set_group) <-  paste0(seqnames(peak_set_group), "-", start(peak_set_group), "-", end(peak_set_group))

# Count number of samples that each peak in group peak set overlap with
sample_presence <- sapply(
  peak_sets,
  function(gr) countOverlaps(peak_set_group, gr) > 0
)

# sample_presence is a logical matrix:
# rows = peaks
# cols = samples

sample_counts <- rowSums(sample_presence)
mcols(peak_set_group)$n_samples <- sample_counts

# Save genomic ranges with all info
save(peak_set_group, file= paste0(subsetDir,cf,"_5e-2.narrowPeak.rdata"))
