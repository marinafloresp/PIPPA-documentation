########################
##       Figure 1     ##
########################

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

## plotting colours/themes ##
col_vals <- c("SOX10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")

pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )
  

cb <- c(
  "#0072B2",  # blue
  "#E69F00",  # orange
  "#F0E442",  # green
  "#CC79A7",  # purple
  "#999999",  # vermillion
  "#56B4E9"   # sky blue
)

pass_fail_cols <- c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
)


## gather data ##
passQC1 <- read.csv(paste0(metaDir,qc1pass),header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)
pheno$fraction[pheno$fraction == "SOX"] <- "Sox10+"
pheno$fraction[pheno$fraction == "NEUN"] <- "NeuN+"
pheno$fraction[pheno$fraction == "IRF8"] <- "IRF8+"

qc1 <- read.csv(file.path(qc1pass), sep = ",",stringsAsFactors = FALSE)
qc1$QCPASS <- as.factor(qc1$QCPASS)
pheno$qcpass1 <- qc1[match(pheno$sampleID,qc1$sampleID),]$QCPASS

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

## plot nucleosome periodicity
propNucleosomesAll.2 <-as.data.frame(propNucleosomesAll)
colnames(propNucleosomesAll.2) <- c("NFR", "Mono", "Di", "Tri", "Quad")
propNucleosomesAll.2$fraction <- pheno[match(rownames(propNucleosomesAll.2), pheno$sampleID),]$fraction
propNucleosomesAll.2$cohort <- pheno[match(rownames(propNucleosomesAll.2), pheno$sampleID),]$cohort
propNucleosomesAll.2$sequencingBatch <- pheno[match(rownames(propNucleosomesAll.2), pheno$sampleID),]$sequencingBatch
propNucleosomesAll.2$passQC1 <- pheno[match(rownames(propNucleosomesAll.2), pheno$sampleID),]$qcpass1
propNucleosomesAll.2$passQC1 <- as.factor(propNucleosomesAll.2$passQC1)

melt.df  <- propNucleosomesAll.2 %>%
pivot_longer(
    cols = 1:5,                     # first 4 columns
    names_to = "FragSize",     # new column for fragment names
    values_to = "percentage"        # new column for values
  )
melt.df$FragSize <- as.factor(melt.df$FragSize)

melt.df$FragSize <- factor(
  melt.df$FragSize,
  levels = c("NFR", "Mono", "Di", "Tri", "Quad")
)

# ########################## Fig 1.a1 ###########################################
# Violin plot for fragment size % (y-axis)
# fragment size group (x-axis)

p1a1 <- ggplot(melt.df,aes(x=FragSize, y=percentage,fill = FragSize))+
   geom_violin(scale = "width",
            position = position_dodge(width = 0.9),
            trim = FALSE) +
  geom_boxplot(width = 0.15,
             position = position_dodge(width = 0.9),
             outlier.shape = NA)+
    labs(y = "Proportion of reads", x = "Fragment size classification", title = "A")+
   stat_summary(
    fun = median,
    geom = "point",
    color = "white",
    size = 2
  )+scale_fill_manual(values = cb)+
      pub_theme +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold")
  )

# ########################## Fig 1.a2 ###########################################
# Violin plot for fragment size % (y-axis)
# fragment size group (x-axis)
# coloured by overall pass/fail Stage 1 QC

p1a2 <- ggplot(melt.df,aes(x=FragSize, y=percentage,fill = passQC1))+
   geom_violin(scale = "width",
            position = position_dodge(width = 0.9),
            trim = FALSE) +
  geom_boxplot(width = 0.15,
             position = position_dodge(width = 0.9),
             outlier.shape = NA)+
    scale_fill_manual(values = pass_fail_cols)+ #c("FALSE" = "#CD5555","TRUE" = "#008B00")) +
    labs(y = "Proportion of reads", x = "Fragment size classification", fill = "Pass Stage 1 QC")+
   stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.9),
    color = "white",
    size = 2
  )+
      pub_theme +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",legend.key.size = unit(0.5, 'cm'),
    plot.title = element_text(size = 16, face = "bold")
  )

## Load Figure 1.b-g made with script Figure1b_cellTypeCheck ##
p1b.g <- readRDS("/fig1bg_ctCheck.rds")

figure1 <- (p1a1 | p1a2) / p1b.g + plot_layout(height=c(1,1.5))

ggsave(
  filename = "Figure_1a.pdf",
  plot = figure1,
  width = 14,
  height = 14
)

ggsave(
  filename = "Figure_1a.png",
  plot = figure1,
  width = 14,
  height = 14
)
