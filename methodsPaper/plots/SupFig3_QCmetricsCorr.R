suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(knitr)
  library(vioplot)
  library(corrplot)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(cowplot)
  library(kableExtra)
  library(scales)
  library(tibble)
  library(ChIPQC)
  library(tidyr)
  library(patchwork)
}))

## for plotting ##
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )
col_vals <- c("SOX10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

## Gather all QC Stage 1 data and results ##
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$fraction[pheno$fraction == "SOX"] <- "SOX10+"
pheno$fraction[pheno$fraction == "NEUN"] <- "NeuN+"
pheno$fraction[pheno$fraction == "IRF8"] <- "IRF8+"
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)

processSum <- read.csv(file.path(paste0(metaDir,"/summariseSampleProcessingProgress.csv")), stringsAsFactors = FALSE, strip.white = TRUE)
fastqc<-read.table(file.path(paste0(fastQCDir,"/multiqc/multiqc_data/multiqc_fastqc.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
readCounts<-read.table(file.path(paste0(alignedDir, "/countMTReads.txt")), fill = TRUE, skip = 1)
alignQC<-read.table(file.path(paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
flagstat<-read.csv(file.path(paste0(alignedDir, "/ENCODEMetrics/collateFlagStatMetrics.txt")), header = FALSE)
flagstat<- flagstat[,c(1:ncol(flagstat)-1)]

missingCol<-setdiff(c("sequencingBatch", "sampleID",  "cohort", "fraction", "experiment", "individualID"), colnames(pheno))
if (length(missingCol) != 0 ) {
  if ("sequencingBatch" %in% missingCol){
    pheno$sequencingBatch<-pheno$cohort
      if (length(missingCol) > 1) {
      warning('Incorrect column names. Missing column(s): ', list(missingCol))
      stop()
    }
  }
}

## exclude duplicates from metadata
dups<-names(which(table(paste(pheno$individualID, pheno$fraction, pheno$tissue, sep = "_")) > 1))
pheno<-unique(pheno[,intersect(c("sequencingBatch", "sampleID","cohort","fraction","individualID", "tissue"), colnames(pheno))])
processSum <- processSum[match(pheno$sampleID, processSum$sampleID),]

## wrangle multiqc stats to retain important columns
colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
mergeStats<-cbind(fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep], fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep])

## calculate percentage of mitochondrial reads
readCounts<-readCounts[match(pheno$sampleID, readCounts$V1),]
mergeStats<-cbind(mergeStats, readCounts$V2/readCounts$V3*100)
colnames(mergeStats)[ncol(mergeStats)]<-"PercentMTReads"

## select alignment columns and match to sampleIDs
aIndex<-match(processSum$sampleID, gsub("\\.bowtie", "", alignQC$Sample))
alignQC$aligned_total<-rowSums(alignQC[,c("paired_aligned_one", "paired_aligned_multi")])
mergeStats<-cbind(mergeStats,alignQC[aIndex,c("overall_alignment_rate", "aligned_total")])
mergeStats<-cbind(sample=processSum$sampleID, mergeStats)

## load flagstat metrics calculated as part of encode qc pipeline
flagstat<-flagstat[match(processSum$sampleID, flagstat$V1),]
files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".pbc.qc")
eMetrics<-NULL
for(each in files){
  tmp<-read.table(paste0(alignedDir, "/ENCODEMetrics/", each))
  if(ncol(tmp) == 7){
    eMetrics<-rbind(eMetrics, tmp)
  } else {
    eMetrics<-rbind(eMetrics, rep(NA, 7))
  }
}
colnames(eMetrics)<-c("TotalReadPairs","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2")
eMetrics$Sample <- files
eMetrics$Sample <- gsub(".pbc.qc", "", eMetrics$Sample)
eIndex<-match(processSum$sampleID, gsub(".pbc.qc", "", files))
eMetrics<-eMetrics[eIndex,]
eMetrics$Sample <- processSum$sampleID
mergeStats<-cbind(mergeStats,eMetrics)

## Option to use a single csv file with all results
fripStats <- read.csv(file.path(paste0(peakDir,"/QCOutput/FRIP_all_samples.csv")), stringsAsFactors = FALSE, strip.white = TRUE)
fripStats$FripMACS3PE <- fripStats$ReadsinBAMPeaks/fripStats$BAMTotalReads
fripStats[is.na(fripStats)] <- 0
fripStats <- fripStats[match(pheno$sampleID, fripStats$SampleName),]

## load periodicity metrics 
propNucleosomesAll<-NULL
diptestStatsAll<-NULL
filePaths <- list.files(file.path(alignedDir,"/QCOutput/"), pattern = "FragmentDistribution_Batch|FragmentDistribution_Batch_")
for(each in filePaths){
  load(paste0(alignedDir,"/QCOutput/", each))
  propNucleosomesAll<-rbind(propNucleosomesAll, propNucleosomes)
  diptestStatsAll<-rbind(diptestStatsAll, diptestStats)
}

propNucleosomesAll<-propNucleosomesAll[match(processSum$sampleID, rownames(propNucleosomesAll)),]

## count number of samples with X million reads
readThres<-seq(0,max(mergeStats[,2], na.rm = TRUE)+10^6, 10^6)
nSamples<-matrix(data = NA, nrow = length(readThres), ncol = length(unique(pheno$sequencingBatch))+1)
for(i in readThres){
  nSamples[1+(i/10^6),1] <- sum(mergeStats[,2] > i, na.rm = TRUE)
  colNum<-2
  for(each in unique(pheno$sequencingBatch)){
    nSamples[1+(i/10^6),colNum] <- sum(mergeStats[which(pheno$sequencingBatch == each),2] > i, na.rm = TRUE)
    colNum<-colNum+1
  }
}

tol <- 0
decreasingProps<-cbind(propNucleosomesAll[,1]+tol > propNucleosomesAll[,2], 
                       propNucleosomesAll[,2]+tol > propNucleosomesAll[,3], 
                       propNucleosomesAll[,3]+tol > propNucleosomesAll[,4]
)
colnames(decreasingProps)<-c("NFR>Mono", "Mono>Di", "Di>Tri")
decreasingProps <- as.data.frame(decreasingProps)

## Load TSS enrichment scores
tssSamples <- read.table(paste0(alignedDir,"/QCOutput/TSS_enrichment.tsv"),header= TRUE)

qcStats <- read.csv(file.path(qc1stats),header = TRUE, sep = ',', stringsAsFactors = FALSE)
qcStats$qcpass <- QCPASS.df[match(qcStats$sampleID, QCPASS.df$sampleID),]$QCPASS
mergeStats$totalreads<-mergeStats$Total.Sequences_R1+mergeStats$Total.Sequences_R2
mergeStats$total_deduplicated_percentage_mean<-(mergeStats$total_deduplicated_percentage_R1+mergeStats$total_deduplicated_percentage_R2)/2
qcStats$NumReadsbeforeQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$totalreads
qcStats$UniqueNumReadsQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$total_deduplicated_percentage_mean

# ########################## Sup Fig 3 ###########################################
# Correlation heatmap of QC metrics refined to those that are used to filter samples and QCPASS indicator variable.

corMergeStats<-cbind(qcStats$UniqueReads,
                     qcStats$NRF,
                     qcStats$PBC1,
                     qcStats$PBC2,
                     qcStats$NFR,
                     qcStats$MNC,
                     qcStats$AlignmentRate)
head(corMergeStats)                
colnames(corMergeStats)<- c('Total\nreads', 'NFR', 'PBC1', 'PBC2', "%NFR", "%Mono", 'Alingment\nrate')

corMergeStats<- corMergeStats[, colSums(corMergeStats, na.rm = TRUE) != 0]

png("SupplementaryFig_3.png", width = 2100, height = 2000, res = 300)

corrplot(
 cor(corMergeStats,use = "pairwise.complete.obs"),
  method = "circle",
  addCoef.col = "black",
  number.cex = 0.8,
  col = colorRampPalette(c("#D55E00", "#F7F7F7", "#0072B2"))(200),
  tl.col = "black",
  tl.cex = 1,
  cl.cex = 1,
  mar = c(1,1,1,1),
  diag = FALSE
)

dev.off()

pdf("SupplementaryFig_3.pdf",width = 11, height = 10)

corrplot(
 cor(corMergeStats, use = "pairwise.complete.obs"),
  method = "circle",
  addCoef.col = "black",
  number.cex = 1.2,
  col = colorRampPalette(c("#D55E00", "#F7F7F7", "#0072B2"))(200),
  tl.col = "black",
  tl.cex = 1,
  cl.cex = 1,
  mar = c(1,1,1,1),
  diag = FALSE
)
dev.off()
