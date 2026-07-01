## ================================================================================================##
##                          ATAC-seq pipeline STEP 7.3: Counts in peaks                            ##
## ================================================================================================##
## EXECUTION: Rscript ./countsInPeaksGroup.r <project>                                             ||
## - execute from scripts directory                                                                ||
##                                                                                                 ||
## DESCRIPTION: This scripts gets read counts in peaks called at cell group level. Only samples    ||
## that passed both QC stages are used, as done to call peaks.                                     ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <configFile> : config file of project on which analysis is being run                          ||
## - peaks called for each cell-type, ending in ".sorted.narrowPeak.filt"                          ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  counts in peaks (all peaks) file: groupCountsAll.rdata in 5_countPeaks folder                ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R version > 4.3, R libraries: Rsubread, GenomicRanges, ChIPpeakAnno,ChIPseeker                ||
##   TxDb.Hsapiens.UCSC.hg38.knownGene,org.Hs.eg.db                                                ||
## - A csv file in the metadata folder with samples used for group peak calling:                   || 
##   samplesGroupAnalysis.csv. This is produced in STEP 7.0                                        ||
## - Parameters in config.r file: groupAnalysisSamples, sampleSheet, macs_thres_group, groupsPeaks ||
##   groupCountsAll,                                                                               ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(GenomicRanges)
  library(ChIPpeakAnno)
  library(ChIPseeker)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library("BSgenome.Hsapiens.UCSC.hg38")
  library(org.Hs.eg.db)  
  library(DESeq2)
  library(ggplot2)
  library(gridExtra)
  library(ggpubr)
  library(edgeR)
}))

## ===============##
##   PARAMETERS   ##
## ===============## 

groupsPeaks <- "/dir/to/peaks/cellfraction/groups"
macs_thres_group <- "threshold used for macs3 peak calling"
groupAnalysisSamples <- "/dir/to/sampleSheet/with/samples/passed/Stage1and2/QC.csv"
groupCountsAll <- "/dir/and/name/for/counts/in/peaks/file.rdata"
alignedDir <-"/dir/to/aligned/bam/files"

## ===============##
##   GATHER DATA  ##
## ===============##

pheno <- read.csv(file.path(groupAnalysisSamples), sep = ",",stringsAsFactors = FALSE)
print(table(pheno$fraction))
print(table(pheno$sequencingBatch))

peakFiles<-list.files(groupsPeaks, pattern=paste0("*_",macs_thres_group,".sorted.narrowPeak.filt"))
print("Peaks files: ")
print(peakFiles)

peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(file.path(groupsPeaks, each)), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
}

peaks<-reduce(unlist(GRangesList(peaks)),min.gapwidth = 150)

## consensus peak set exported
export(peaks, con = paste0(groupsPeaks,"/consensus4CT_",macs_thres_group,".narrowPeak.bed"), format = "bed")

## Annotate peaks
peaksAnno<-annotatePeak(peaks, tssRegion=c(-1000, 1000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", verbose = FALSE)

bamFiles<-list.files(alignedDir, pattern = "*.filt.nodup.bam$", recursive = TRUE, include.dirs = TRUE)
sampleNames<-gsub(".filt.nodup.bam$", "", bamFiles)
sampleNames <- sampleNames[sampleNames %in% pheno$sampleID]
bamFiles<-bamFiles[bamFiles %in% unlist(lapply(sampleNames, function(x) paste0(x,".filt.nodup.bam")))]
print(bamFiles)
                                               
print(dim(bamFiles))
stopifnot(all(dim(bamFiles) == nrow(pheno)))

## ================##
##   COUNTS (ALL)  ##
## ================##

## File name specified in config file
print(groupCountsAll)

ctPeaks <- as.data.frame(peaksAnno@anno)
ctPeaks$GeneID <- paste0("Peak_", row.names(ctPeaks))
ctPeaks <- ctPeaks[,c("GeneID", "seqnames", "start", "end", "strand")]
colnames(ctPeaks) <- c("GeneID", "Chr", "Start", "End", "Strand")
## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = groupCountsAll)
