      ########################################################
## Script to plot fragment size distribution for pass/fail samples ##
      ########################################################
args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)
library(ggplot2)
library(ATACseqQC)
library(plyr)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(reshape2)
library(tidyr)

## parameters ##
outputDir <- "/dir/to/plot/file"
qc1pass <- "/dir/to/qc/stage/1/output/file.csv"
sampleSheet <- "/dir/to/sample/sheet.csv"
alignedDir <- "/dir/to/aligned/dir"

## gather data ##
passQC1 <- read.csv(paste0(metaDir,qc1pass),header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)
pheno$fraction[pheno$fraction == "SOX"] <- "Sox10+"
pheno$fraction[pheno$fraction == "NEUN"] <- "NeuN+"
pheno$fraction[pheno$fraction == "IRF8"] <- "IRF8+"

## gather fragment size distribution files produced at STEP 2 of QC pipeline ##
processSum <- pheno
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
diptestStatsAll <- as.data.frame(diptestStatsAll)
diptestStatsAll$sampleCode <- processSum$sampleCode
propNucleosomesAll <-as.data.frame(propNucleosomesAll)
propNucleosomesAll$fraction <- pheno[match(rownames(propNucleosomesAll), pheno$sampleID),]$fraction

## Plot ##
setwd(outputDir)
pdf(paste0("fragSize.pdf"), width = 10, height = 10)
head(propNucleosomesAll)

## For only NEUN,SOX10, split by pass/fail
propNucleosomesAll.1<- propNucleosomesAll[propNucleosomesAll$fraction %in% c("NeuN+","Sox10+"),]
propNucleosomesAll.1$passQC1 <- as.factor(propNucleosomesAll.1$passQC1)
head(propNucleosomesAll.1)
melt.df  <- propNucleosomesAll.1 %>%
pivot_longer(
    cols = 1:5,                     # first 4 columns
    names_to = "FragSize",     # new column for fragment names
    values_to = "percentage"        # new column for values
  )
head(melt.df)

ggplot(melt.df,aes(x=as.factor(FragSize), y=percentage))+
   geom_violin(aes(fill = passQC1),scale = "width",
            position = position_dodge(width = 0.9),
            trim = FALSE) +
  geom_boxplot(width = 0.15,
             position = position_dodge(width = 0.9),
             outlier.shape = NA)+
    facet_grid(fraction ~ passQC1)+
    scale_fill_manual(values = c("FALSE" = "#CD5555","TRUE" = "#008B00")) +
    labs(y = "Proportion of reads", x = "Fragment size classification", fill = "Pass Periodicity (Stage 1)")+
    stat_summary(fun="median", geom="point", color="white")+
    scale_x_discrete(labels = c("NFR", "Mono", "Di", "Tri", "Quad"))+
    theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  )

## For all cell-types, split only by pass/fail
propNucleosomesAll.2<- propNucleosomesAll[propNucleosomesAll$fraction %in% c("NeuN+","Sox10+"),]
propNucleosomesAll.2$passQC1 <- as.factor(propNucleosomesAll.2$passQC1)
melt.df  <- propNucleosomesAll.2 %>%
pivot_longer(
    cols = 1:5,                     # first 4 columns
    names_to = "FragSize",     # new column for fragment names
    values_to = "percentage"        # new column for values
  )
ggplot(melt.df,aes(x=as.factor(FragSize), y=percentage))+
   geom_violin(aes(fill = passQC1),scale = "width",
            position = position_dodge(width = 0.9),
            trim = FALSE) +
  geom_boxplot(width = 0.15,
             position = position_dodge(width = 0.9),
             outlier.shape = NA)+
    facet_grid(rows = vars(passQC1))+
    scale_fill_manual(values = c("FALSE" = "#CD5555","TRUE" = "#008B00")) +
    labs(y = "Proportion of reads", x = "Fragment size classification", fill = "Pass Periodicity (Stage 1)")+
    stat_summary(fun="median", geom="point", color="white")+
    scale_x_discrete(labels = c("NFR", "Mono", "Di", "Tri", "Quad"))+
    theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  )
dev.off()
propNucleosomesAll$passQC1 <- passQC1[match(rownames(propNucleosomesAll), passQC1$sampleID),]$Periodicity
