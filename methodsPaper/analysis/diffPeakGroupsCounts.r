## =============================================================================================================##
##                         ATAC-seq pipeline STEP 7.4: Differential analysis groups peaks                       ##
## =============================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/diffPeakGroupsCounts.r                                      ||
## - execute from scripts directory                                                                             ||
##                                                                                                              ||
## DESCRIPTION: Performs differential accessible regions analysis between cell-groups based on counts in peaks  //
##  in order to perform cell-type check.                                                                        ||
##                                                                                                              ||
## INPUTS:                                                                                                      ||
## - <project> : project on which analysis is being run                                                         ||
##                                                                                                              ||
## OUTPUTS:                                                                                                     ||
## - cellType_diffPeaks_edgeR_macs_thres_group.filt.narrowPeak.rdata: differential analysis results in folder   ||
##   for group counts.                                                                                          ||
##                                                                                                              ||
## REQUIRES:                                                                                                    ||
## - R/4.2.1-foss-2022a, libraries: Rsubread, GenomicRanges, ChIPpeakAnno, ChIPseeker ,ggplot2, ggpubr          ||
## TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db                                                              ||
## - Rdata file with counts in peaks results from step 7.3                                                      ||
## - Parameters specified in config.r file: "groupAnalysisSamples","sampleSheet", "groupsCounts", "refCellType",||
##   "groupCountsAll"                                                                                           || 
## =============================================================================================================##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(GenomicRanges)
  library(RColorBrewer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(edgeR)
  library(ggplot2)
  library(ggpubr)
  library(DESeq2)
}))

## Load counts in peaks files and match to pheno samples data
loadCounts <- function(counts){
  colnames(counts) <- lapply(colnames(counts),function(x) gsub("*.filt.nodup.bam", "",x))
  counts <- counts[, colSums(is.na(counts)) == 0]
  return(counts)
}

## ===============##
##   PARAMETERS   ##
## ===============## 

groupsPeaks <- "/dir/to/peaks/cellfraction/groups"
macs_thres_group <- "threshold used for macs3 peak calling"
groupAnalysisSamples <- "/dir/to/sampleSheet/with/samples/passed/Stage1and2/QC.csv"
groupCountsAll <- "/dir/and/name/for/counts/in/peaks/file.rdata"
alignedDir <-"/dir/to/aligned/bam/files"
diffPeaksCountsAll <-"/dir/and/name/for/diff/counts/in/peaks/file.rdata"
refCellType <-"cell-type to use as reference"
                             
## ===============##
##   GATHER DATA  ##
## ===============##

pheno <- read.csv(file.path(groupAnalysisSamples), sep = ",",stringsAsFactors = FALSE)
                             
## Load counts in peaks file
load(groupCountsAll)
counts <- loadCounts(fc_ctPeaks$count)
cnt_table <- fc_ctPeaks$annotation
pheno.meta <- pheno[pheno$sampleID %in% colnames(counts),]
counts <- counts[,match(pheno.meta$sampleID,colnames(counts))]
dim(counts)
dim(pheno.meta)
                             
print("Number of samples per cell-fraction")
table(pheno.meta$fraction)

dge <- DGEList(counts=counts,
           group=pheno.meta$fraction)
keep <- filterByExpr(dge$counts, group = pheno.meta$fraction)
dge <- dge[keep, , keep.lib.sizes = FALSE]

pheno.meta$fraction <- as.factor(pheno.meta$fraction)
                             
## refCellType: base level cell-fraction should be specified in config file
pheno.meta$fraction = relevel(pheno.meta$fraction, ref = refCellType)
design <- model.matrix(~ fraction , data =pheno.meta )
d <- DGEList(counts=dge$counts,
           group=pheno.meta$fraction)
d <- calcNormFactors(d)
d <- estimateDisp(d, design, robust = TRUE)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2) 
results <- lrt$table
results <- as.data.frame(results[order(results$PValue),])
results$Geneid <- rownames(results)
results$Bon <- p.adjust(results$PValue, method="bonferroni")

## File name specified in config file
save(results,  file = diffPeaksCountsAll)




