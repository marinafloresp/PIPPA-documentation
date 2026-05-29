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


args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)
col_vals <- c("Sox10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")

## Gather all QC Stage 1 data and results ##
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$fraction[pheno$fraction == "SOX"] <- "Sox10+"
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

#diptestStatsAll<-as.data.frame(diptestStatsAll)
propNucleosomesAll<-propNucleosomesAll[match(processSum$sampleID, rownames(propNucleosomesAll)),]
#diptestStatsAll<-diptestStatsAll[match(processSum$sampleID, sub("X", "", gsub("\\.D$", "", rownames(diptestStatsAll)))),]
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

#propNucleosomesAll <- na.omit(propNucleosomesAll)
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
tssSamples <- read.table(paste0(alignedDir,"/QCOutput/TSS_enrichment.tsv"),header= TRUE)

# ########################## Sup Fig 2.a ###########################################
# Histogram of number of reads (pre QC) 

mergeStats <- mergeStats[,-c(1,17)]

dfa <- data.frame(pheno$sampleID, mergeStats[,c(1,9)]) %>%
  setNames(c('Samples', 'Raw reads')) %>%
  melt(id='Samples') %>% as.data.frame()

plot.a <- ggplot(dfa, aes(x=Samples, y=value)) +
  geom_bar(position='dodge', stat='identity') +
  labs(x = "", y='Number of reads (x10^6)', title="A") +
  geom_hline(color='red', size=0.8, aes(yintercept=2)) +
  pub_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position='bottom',
    legend.title=element_blank()
  ) +
  scale_fill_manual(values=colorBlindGrey8) +
  scale_y_continuous(
    labels = label_number(suffix="M", scale=1e-6)
  )

# ########################## Sup Fig 2.b ###########################################
# plot %unique reads against read number

plot.b <- ggplot(
  mergeStats,
  aes(x = Sample, y = total_deduplicated_percentage)
)+ geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = '% unique reads', title = "B")+
  geom_hline(color='red', size=0.8, aes(yintercept=2)) +
  pub_theme + 
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom')

# ########################## Sup Fig 2.c ###########################################
# Histogram of alignment rate 

df.c <- mergeStats[,c('aligned_total','overall_alignment_rate',"Sample")]
plot.c <- ggplot(df.c, aes(x = Sample, y = overall_alignment_rate))+
  geom_bar(position = 'dodge', stat = 'identity')+
  geom_hline(yintercept=80, #acceptable
             color ='red')+
  labs(x = '', y = '% Alignment rate', title = "C")+
  pub_theme + 
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom')

# ########################## Sup Fig 2.d ###########################################
# Histogram of number of reads (post QC) 

df.d <- mergeStats[,c('aligned_total','overall_alignment_rate',"Sample")]
plot.d <- ggplot(df.d, aes(x = Sample, y = aligned_total))+
  geom_bar(position = 'dodge', stat = 'identity')+
  geom_hline(yintercept=c(95, 80), #acceptable
             color = c('black', colorBlindGrey8[2]))+
  labs(x = '', y = 'Number of reads aligned', title = "D")+
  pub_theme + 
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom')+
	scale_y_continuous(labels = label_number(suffix = "M", scale = 1e-6))

# ########################## Sup Fig 2.e ###########################################
# Histogram of % MT reads

df.d <- mergeStats[,c('PercentMTReads',"Sample")]
plot.e <- ggplot(df.d, aes(x = Sample, y = PercentMTReads))+
  geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = '% MT reads', title = "E")+
  pub_theme + 
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')

# ########################## Sup Fig 2.f ###########################################
# Histogram of ENDCODE metrics (PBC1, PBC2, NRF) with thresholds labelled 

plot.f <- ggplot(eMetrics, aes(x = Sample, y = PBC1))+
  geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = 'PBC1')+ 
  pub_theme + geom_hline(yintercept=pbc1Thres, #acceptable
             color = "red")+
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')
        
plot.g <- ggplot(eMetrics, aes(x = Sample, y = PBC2))+
  geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = 'PBC2')+ 
  pub_theme + geom_hline(yintercept=pbc2Thres, #acceptable
             color = "red")+
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')        

plot.h <- ggplot(eMetrics, aes(x = Sample, y = NRF))+
  geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = 'NRF')+ 
  pub_theme + geom_hline(yintercept=nrfThres, #acceptable
             color = "red")+
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')   

# ########################## Sup Fig 2.g ###########################################
# Histogram of TSS enrichment distribution

plot.i<- ggplot(tssSamples, aes(x = sample, y =TSS_enrichment)) + 
  geom_bar(position = 'dodge', stat = 'identity')+
  labs(x = '', y = 'TSS enrichment score', title = "G")+ 
  pub_theme + 
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')

########################### Sup Fig 2 ######################################
### Gather all plots 

design <- "AAABBB
           AAABBB
           CCCDDD
           CCCDDD
           EEEIII
           EEEIII
           FFGGHH
           FFGGHH"

combined_fig2 <- plot.a + plot.b + plot.c + 
  plot.d + plot.e + plot.f + plot.g+ plot.h + plot.i + 
  plot_layout(design = design)

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
