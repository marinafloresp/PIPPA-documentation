######################################
##       Supplementary Fig 2,3      ##
######################################

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
  library(plotly)
  library(tidyr)
  library(pheatmap)
  library(cowplot)
  library(patchwork)
}))

## load config variables
args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
cf<-args[2]
source(configFile)

col_vals_pass <- c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
)

## Gather all QC Stage 1 data and results ##
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)

processSum <- read.csv(file.path(paste0(metaDir,"/summariseSampleProcessingProgress.csv")), stringsAsFactors = FALSE, strip.white = TRUE)
fastqc<-read.table(file.path(paste0(fastQCDir,"/multiqc/multiqc_data/multiqc_fastqc.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
readCounts<-read.table(file.path(paste0(alignedDir, "/countMTReads.txt")), fill = TRUE, skip = 1)
alignQC<-read.table(file.path(paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
flagstat<-read.csv(file.path(paste0(alignedDir, "/ENCODEMetrics/collateFlagStatMetrics.txt")), header = FALSE)
flagstat<- flagstat[,c(1:ncol(flagstat)-1)]

# wrangle data
missingCol<-setdiff(c("sequencingBatch", "sampleID", "sampleCode", "cohort", "fraction", "experiment", "individualID"), colnames(pheno))
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
pheno<-unique(pheno[,intersect(c("sequencingBatch", "sampleID","sampleCode", "cohort","fraction","individualID", "tissue"), colnames(pheno))])

idx <- match(pheno$sampleID, processSum$sampleID)
processSum <- processSum[idx[!is.na(idx)], ]
processSum$sampleCode <- pheno[match(processSum$sampleID, pheno$sampleID),]$sampleCode
colsKeep<-c("sampleID","sampleCode", "R1Filename","R2Filename")
mergeStats <- processSum[,colsKeep]

## wrangle multiqc stats to retain important columns
colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
R1stats <- fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep]
colnames(R1stats) <- paste0(colnames(R1stats) , "_R1")
R2stats <- fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep]
colnames(R2stats) <- paste0(colnames(R2stats) , "_R2")
mergeStats<-cbind(mergeStats,R1stats,R2stats )

## calculate percentage of mitochondrial reads
readCounts<-readCounts[match(mergeStats$sampleID, readCounts$V1),]
MTpercentReads <- readCounts$V2/readCounts$V3*100
mergeStats$"PercentMTReads" <- MTpercentReads

## select alignment columns and match to sampleIDs
aIndex<-match(mergeStats$sampleID, gsub("\\.bowtie", "", alignQC$Sample))
alignQC$aligned_total<-rowSums(alignQC[,c("paired_aligned_one", "paired_aligned_multi")])
mergeStats<-cbind(mergeStats,alignQC[aIndex,c("overall_alignment_rate", "aligned_total")])

## load flagstat metrics calculated as part of encode qc pipeline
files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = "\\.pbc\\.qc$")
eMetrics <- do.call(rbind, lapply(files, function(each) {
  tmp <- read.table(
    file.path(alignedDir, "ENCODEMetrics", each)
  )
  if (ncol(tmp) != 7) {
    tmp <- as.data.frame(matrix(NA, nrow = 1, ncol = 7))
  }
  tmp$sampleID <- sub("\\.pbc\\.qc$", "", each)
  tmp
}))

colnames(eMetrics)<-c("TotalReadPairs_afterAlignment","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2", "Sample")
eIndex<-match(mergeStats$sampleID, eMetrics$Sample)
eMetrics<-eMetrics[eIndex,]
mergeStats<-cbind(mergeStats,eMetrics)
mergeStats <- subset(mergeStats, select = -Sample)

## Merge Peak calling and FRIP results
fripStats <- read.csv(file.path(fripFile), stringsAsFactors = FALSE, strip.white = TRUE)
fripStats$FripMACS3PE <- fripStats$ReadsinBAMPeaks/fripStats$BAMTotalReads
eIndex<-match(mergeStats$sampleID, fripStats$SampleName)
fripStats<-fripStats[eIndex,]
mergeStats <- cbind(mergeStats, fripStats[, c("ReadsinBAMPeaks","BAMTotalReads","BAMPeaks","FripMACS3PE")] )
mergeStats$TotalNumReads <- mergeStats$BAMTotalReads

## load periodicity metrics 
filePaths <- list.files(file.path(alignedDir,"/QCOutput/"), pattern = "FragmentDistribution_*")
propNucleosomesAll <- do.call(rbind, lapply(filePaths, function(each) {

  load(file.path(alignedDir, "QCOutput", each))

  propNucleosomes
}))
colnames(propNucleosomesAll) <- c("NFR", "Mono", "Di", "Tri", "Quad")
eIndex<-match(mergeStats$sampleID, rownames(propNucleosomesAll))
propNucleosomesAll<-propNucleosomesAll[eIndex,]
mergeStats <- cbind(mergeStats,propNucleosomesAll)
                    
diptestStatsAll <- do.call(rbind, lapply(filePaths, function(each) {

  load(file.path(alignedDir, "QCOutput", each))

  diptestStats
}))
eIndex<-match(mergeStats$sampleID, gsub("\\.D$", "", rownames(diptestStatsAll)))
diptestStatsAll<-diptestStatsAll[eIndex,]
mergeStats <- cbind(mergeStats,diptestStatsAll)

## add more metadata info
mergeStats <- mergeStats %>%
  left_join(pheno %>% select("sampleID","sequencingBatch","cohort", "fraction"), by = "sampleID")

## count number of samples with X million reads
readThres<-seq(0,max(mergeStats$Total.Sequences_R1, na.rm = TRUE)+10^6, 10^6)
nSamples<-matrix(data = NA, nrow = length(readThres), ncol = length(unique(mergeStats$sequencingBatch))+1)
for(i in readThres){
  nSamples[1+(i/10^6),1] <- sum(mergeStats$Total.Sequences_R1 > i, na.rm = TRUE)
  colNum<-2
  for(each in unique(mergeStats$sequencingBatch)){
    nSamples[1+(i/10^6),colNum] <- sum(mergeStats[which(mergeStats$sequencingBatch == each),3] > i, na.rm = TRUE)
    colNum<-colNum+1
  }
}

## a successful ATAC experiment should have decreasing proportions in the nucleosome, mono, di,tri- nucleosomes
## define a tolerance for these comparisons
tol <- 0
decreasingProps<-cbind(mergeStats$NFR+tol > mergeStats$Mono, 
                       mergeStats$Mono+tol > mergeStats$Di, 
                       mergeStats$Di+tol > mergeStats$Tri
)
colnames(decreasingProps)<-c("NFR>Mono", "Mono>Di", "Di>Tri")
mergeStats <- cbind(mergeStats,decreasingProps)

## Load TSS enrichment scores
tssSamples <- read.table(tssEnrich,header= TRUE)
eIndex<-match(mergeStats$sampleID, tssSamples$sample)
tssSamples<-tssSamples[eIndex,]
mergeStats <- cbind(mergeStats,TSSScore = tssSamples$"TSS_enrichment")

QCPASS.df <- read.csv(file.path(qc1pass),header = TRUE, sep = ',', stringsAsFactors = FALSE)
QCPASS.df <- QCPASS.df[,-1]
qcStats <- read.csv(file.path(qc1stats),header = TRUE, sep = ',', stringsAsFactors = FALSE)
qcStats$qcpass <- QCPASS.df[match(qcStats$sampleID, QCPASS.df$sampleID),]$QCPASS
mergeStats$totalreads<-mergeStats$Total.Sequences_R1+mergeStats$Total.Sequences_R2
mergeStats$total_deduplicated_percentage_mean<-(mergeStats$total_deduplicated_percentage_R1+mergeStats$total_deduplicated_percentage_R2)/2
qcStats$NumReadsbeforeQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$totalreads
qcStats$UniqueNumReadsQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$total_deduplicated_percentage_mean

head(qcStats)
table(qcStats$QCPASS)

pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.key.width = unit(2.2, "cm")
  )

dfa <- data.frame(
  Samples = qcStats$sampleID,
  qcpass = qcStats$qcpass,
  reads = qcStats$NumReadsbeforeQC
) 
 
plot.a <- ggplot(dfa, aes( x = reads)) +
  geom_histogram(bins = 100,fill ='gray')+
  labs( x = 'Number of reads (R1 + R2)', y = "Number of samples", title = "A") +
  pub_theme +
  scale_x_continuous(
    labels = label_number(suffix = "M", scale = 1e-6)
  )

dfb <- data.frame(
  Samples = qcStats$sampleID,
  reads = qcStats$UniqueNumReadsQC
)
plot.b <- ggplot(dfb,aes(x =  reads))+
geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% unique reads (mean R1 + R2)', y = "Number of samples", title = "B")+
  pub_theme 

align_f <- qcStats
align_f$qcpass<- ifelse(align_f$AlignmentRate > 80, "TRUE", "FALSE")

plot.c <- ggplot(align_f, aes(x = AlignmentRate,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=80, 
             color ='red')+
             scale_fill_manual(values = col_vals_pass) +
  labs(x  = '% Alignment rate', y = "Number of samples", title = "C", fill = "PASS")+
  pub_theme 

mergeStats$qcpass <- ifelse(mergeStats$BAMTotalReads > numReads*10^6, "TRUE", "FALSE")
plot.d <- ggplot(mergeStats, aes(x = BAMTotalReads,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=numReads*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x =  'Number of reads aligned (after QC)', y = "Number of samples", fill = "PASS",title = "D")+
  pub_theme +
	scale_x_continuous(labels = label_number(suffix = "M", scale = 1e-6))

plot.e <- ggplot(mergeStats, aes(x =  PercentMTReads))+
  geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% MT reads', y = "Number of samples", title = "E")+
  pub_theme 
  
eMetrics_f <- qcStats
eMetrics_f$qcpass <- ifelse(eMetrics_f$PBC1 > pbc1Thres, "TRUE", "FALSE")

eMetrics_g <- qcStats
eMetrics_g$qcpass <- ifelse(eMetrics_g$PBC2 > pbc2Thres, "TRUE", "FALSE")

eMetrics_h <- qcStats
eMetrics_h$qcpass <- ifelse(eMetrics_h$NRF > nrfThres, "TRUE", "FALSE")

plot.f <- ggplot(eMetrics_f, aes(x = PBC1, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC1', y = "Number of samples", title = "G", fill = "PASS") +
  pub_theme + geom_vline(xintercept = pbc1Thres, color = "red")

plot.g <- ggplot(eMetrics_g, aes(x = PBC2, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC2', y = "Number of samples", title = "H", fill = "PASS") +
  pub_theme + geom_vline(xintercept = pbc2Thres, color = "red")

plot.h <- ggplot(eMetrics_h, aes(x = NRF, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NRF', y = "Number of samples", title = "I", fill = "PASS") +
  pub_theme + geom_vline(xintercept = nrfThres, color = "red")+ theme(
  legend.position = "none"
)

plot.i<- ggplot(qcStats, aes(x  =TSS_score)) + 
  geom_histogram(bins = 100,fill ='gray')+
  labs(x =  'TSS enrichment score', y = "Number of samples", title = "F")+ 
  pub_theme 

propNuc_j <-qcStats
propNuc_j$qcpass <- ifelse(propNuc_j$NFR > nucfThres, "TRUE", "FALSE")

propNuc_k <-qcStats
propNuc_k$qcpass <- ifelse(propNuc_k$MNC > monoThres, "TRUE", "FALSE")

plot.j <- ggplot(propNuc_j, aes(x = NFR, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NFR fragments', y = "Number of samples", title = "J", fill = "PASS") +
  pub_theme + geom_vline(xintercept = nucfThres, color = "red") + 
  theme(legend.position = "none")

plot.k <- ggplot(propNuc_k, aes(x = MNC, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'MNC fragments', y = "Number of samples", title = "K", fill = "PASS") +
  pub_theme + geom_vline(xintercept = monoThres, color = "red")

plot.c <- plot.c + theme(legend.position = "none")
plot.f <- plot.f + theme(legend.position = "none")

plot.d <- plot.d + theme(legend.position = "none")
plot.g <- plot.g + theme(legend.position = "none")
plot.h <- plot.h + theme(legend.position = "none")
plot.i <- plot.i + theme(legend.position = "none")
plot.j <- plot.j + theme(legend.position = "none")
plot.k <- plot.k + theme(legend.position = "none")


combined_figall <- (plot.a + plot.b + plot.c +
  plot.d + plot.e + plot.f +
  plot.g + plot.h + plot.i + 
  plot.j + plot.k + 
  plot_layout( design = "
AAABBB
AAABBB
CCCDDD
CCCDDD
EEEIII
EEEIII
FFGGHH
FFGGHH
JJJKKK
JJJKKK
"
)+ 
  plot_layout(guides = "collect"             # <-- collect all legends into guide_area()
  ) & theme(
  legend.position = "bottom"
))

ggsave(
  filename = paste0("SupplementaryFig_2.png"),
  plot = combined_figall,
  width = 12,
  height = 17)
  
ggsave(
  filename = paste0("SupplementaryFig_2.pdf"),
  plot = combined_figall,
  width = 12,
  height = 17
)

