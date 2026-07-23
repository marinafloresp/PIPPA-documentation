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
col_vals_pass <- c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
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
diptestStatsAll<-diptestStatsAll[match(processSum$sampleID, gsub("\\.D$", "", rownames(diptestStatsAll))),]
rownames(diptestStatsAll)<-processSum$sampleID

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

## a successful ATAC experiment should have decreasing proportions in the nucleosome, mono, di,tri- nucleosomes
## define a tolerance for these comparisons
tol <- 0
decreasingProps<-cbind(propNucleosomesAll[,1]+tol > propNucleosomesAll[,2], 
                       propNucleosomesAll[,2]+tol > propNucleosomesAll[,3], 
                       propNucleosomesAll[,3]+tol > propNucleosomesAll[,4]
)
colnames(decreasingProps)<-c("NFR>Mono", "Mono>Di", "Di>Tri")
decreasingProps <- as.data.frame(decreasingProps)

## Load TSS enrichment scores
tssSamples <- read.table(tssEnrich,header= TRUE)
eIndex<-match(mergeStats$sampleID, tssSamples$sample)
tssSamples<-tssSamples[eIndex,]
mergeStats <- cbind(mergeStats,TSSScore = tssSamples$"TSS_enrichment")

QCPASS.df <- read.csv(file.path(qc1pass),header = TRUE, sep = ',', stringsAsFactors = FALSE)
QCPASS.df <- QCPASS.df[,-1]
head(QCPASS.df)
qcStats <- read.csv(file.path(qc1stats),header = TRUE, sep = ',', stringsAsFactors = FALSE)
qcStats$qcpass <- QCPASS.df[match(qcStats$sampleID, QCPASS.df$sampleID),]$QCPASS
mergeStats$totalreads<-mergeStats$Total.Sequences_R1+mergeStats$Total.Sequences_R2
mergeStats$total_deduplicated_percentage_mean<-(mergeStats$total_deduplicated_percentage_R1+mergeStats$total_deduplicated_percentage_R2)/2
qcStats$NumReadsbeforeQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$totalreads
qcStats$UniqueNumReadsQC <- mergeStats[match(qcStats$sampleID,mergeStats$sampleID),]$total_deduplicated_percentage_mean

# ########################## Sup Fig 2.a ###########################################
# Histogram of number of reads (pre QC) 

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

# ########################## Sup Fig 2.b ###########################################
# plot %unique reads against read number

dfb <- data.frame(
  Samples = qcStats$sampleID,
  reads = qcStats$UniqueNumReadsQC
)
plot.b <- ggplot(dfb,aes(x =  reads))+
geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% unique reads (mean R1 + R2)', y = "Number of samples", title = "B")+
  pub_theme 

# ########################## Sup Fig 2.c ###########################################
# Histogram of alignment rate 

dfc <- qcStats
dfc$qcpass<- ifelse(dfc$AlignmentRate > 80, "TRUE", "FALSE")

plot.c <- ggplot(dfc, aes(x = AlignmentRate,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=80, 
             color ='red')+
             scale_fill_manual(values = col_vals_pass) +
  labs(x  = '% Alignment rate', y = "Number of samples", title = "C", fill = "PASS")+
  pub_theme 

# ########################## Sup Fig 2.d ###########################################
# Histogram of number of reads (post QC) 

mergeStats$qcpass <- ifelse(mergeStats$BAMTotalReads > uniqReads*2*10^6, "TRUE", "FALSE")
plot.d <- ggplot(mergeStats, aes(x = BAMTotalReads,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=uniqReads*2*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x =  'Number of reads aligned (after QC)', y = "Number of samples", fill = "PASS",title = "D")+
  pub_theme +
	scale_x_continuous(labels = label_number(suffix = "M", scale = 1e-6))

# ########################## Sup Fig 2.e ###########################################
# Histogram of % MT reads

plot.e <- ggplot(mergeStats, aes(x =  PercentMTReads))+
  geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% MT reads', y = "Number of samples", title = "E")+
  pub_theme 

# ########################## Sup Fig 2.f-h ###########################################
# Histogram of ENDCODE metrics (PBC1, PBC2, NRF) with thresholds labelled 

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

# ########################## Sup Fig 2.i ###########################################
# Histogram of TSS enrichment distribution

plot.i<- ggplot(qcStats, aes(x  =TSS_score)) + 
  geom_histogram(bins = 100,fill ='gray')+
  labs(x =  'TSS enrichment score', y = "Number of samples", title = "F")+ 
  pub_theme 

# ########################## Sup Fig 2.j-k ###########################################
# Histogram of %NFR and %MNC

propNuc_j <-qcStats
propNuc_j$qcpass <- ifelse(propNuc_j$NFR > nucfThres, "TRUE", "FALSE")

propNuc_k <-qcStats
propNuc_k$qcpass <- ifelse(propNuc_k$MNC > monoThres, "TRUE", "FALSE")

plot.j <- ggplot(propNuc_j, aes(x = NFR, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = '% NFR fragments', y = "Number of samples", title = "J", fill = "PASS") +
  pub_theme + geom_vline(xintercept = nucfThres, color = "red") + 
  theme(legend.position = "none")

plot.k <- ggplot(propNuc_k, aes(x = MNC, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = '% MNC fragments', y = "Number of samples", title = "K", fill = "PASS") +
  pub_theme + geom_vline(xintercept = monoThres, color = "red")

########################### Sup Fig 2 ######################################
### Gather all plots 

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
  plot_layout(guides = "collect"             
  ) & theme(
  legend.position = "bottom"
))


ggsave(
  filename = "SupplementaryFig_2.png",
  plot = combined_fig2,
  width = 12,
  height = 12
)

ggsave(
  filename = "SupplementaryFig_2.pdf",
  plot = combined_fig2,
  width = 12,
  height = 15
)
