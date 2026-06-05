## Supplementary Figure 4: Stage 2 QC results ## 

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

pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.key.width = unit(2.2, "cm")
  )

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

col_vals <- c("Sox10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")

col_vals_pass <-c("FALSE" = "#E15759", "TRUE" = "#59A14F")

## Read data ## 
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$fraction[pheno$fraction == "SOX"] <- "Sox10+"
pheno$fraction[pheno$fraction == "NEUN"] <- "NeuN+"
pheno$fraction[pheno$fraction == "IRF8"] <- "IRF8+"
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)

qcPass1 <- read.csv(qcStatsStage1.csv, stringsAsFactors = FALSE, strip.white = TRUE, header = TRUE)
qcPass1$QCPASS <- as.factor(qcPass1$QCPASS)

qcPass2 <- read.csv(qcStatsStage2.csv, stringsAsFactors = FALSE, strip.white = TRUE, header = TRUE)
qcPass2$QCPASS <- as.factor(qcPass2$QCPASS)
pheno$qcpass1 <- qcPass1[match(pheno$sampleID,qcPass1$sampleID),]$QCPASS
pheno$qcpass2 <- qcPass2[match(pheno$sampleID,qcPass2$sampleID),]$QCPASS

### Collate sex check results ###
xPeakFile<-file.path(sexChrCountsDir, "/chrX.peakcounts.txt")
yPeakFile<-file.path(sexChrCountsDir, "/chrY.peakcounts.txt")

xPeaks<-read.table(xPeakFile, stringsAsFactors = FALSE)
yPeaks<-read.table(yPeakFile, stringsAsFactors = FALSE,fill = TRUE)

colnames(xPeaks)[c(4,18,19)] <- c("sex-gene","counts","sampleID")
colnames(yPeaks)[c(4,11,12)] <- c("peak-name","counts","sampleID")

xPeaks$sampleID<-basename(xPeaks$"sampleID")
yPeaks$sampleID<-basename(yPeaks$"sampleID")

xPeaksWide = xPeaks %>% 
	select("sex-gene","counts","sampleID") %>% 
  spread("sex-gene","counts")

yPeaksWide = yPeaks %>% 
	select("peak-name","counts","sampleID") %>%
 spread("peak-name","counts")

xPeaksWide$sampleID<-gsub(".chrX.tn5.tagAlign.gz", "", xPeaksWide$sampleID, fixed = TRUE)
yPeaksWide$sampleID<-gsub(".chrY.tn5.tagAlign.gz", "", yPeaksWide$sampleID,fixed = TRUE)

xPeaksWide<-xPeaksWide[match(pheno$sampleID,xPeaksWide$sampleID), ]
yPeaksWide<-yPeaksWide[match(pheno$sampleID,yPeaksWide$sampleID),]

xPeaksWide$YTot<-rowSums(yPeaksWide[,-1])
xPeaksWide$sex<-pheno$gender
xPeaksWide$sampleCode <- pheno[match(xPeaksWide$sampleID,pheno$sampleID),]$sampleCode
xPeaksWide$qcpass2 <- pheno[match(xPeaksWide$sampleID,pheno$sampleID),]$qcpass2


## determine threshold as a line to bisect the two groups 
ratioF<-xPeaksWide$YTot/xPeaksWide$FIRRE
ratioF[!is.finite(ratioF)]<-NA
sexMeans<-aggregate(ratioF, by = list(xPeaksWide$sex), mean, na.rm = TRUE)[,2]
thresF<-min(sexMeans)+ 0.5*abs(diff(sexMeans))

ratioX<-xPeaksWide$YTot/xPeaksWide$XIST
ratioX[!is.finite(ratioX)]<-NA
sexMeans<-aggregate(ratioX, by = list(xPeaksWide$sex), mean, na.rm = TRUE)[,2]
thresX<-min(sexMeans)+ 0.5*abs(diff(sexMeans))

sexPredict<-data.frame("sampleID" = xPeaksWide$sampleID,"LabelledSex" = xPeaksWide$sex, "FIRREratio" = ratioF > thresF,"XISTratio" = ratioX > thresX)

sexPredict[sexPredict == TRUE]<-"M"
sexPredict[sexPredict == FALSE]<-"F"
sexPredict$concordantPredictions<-sexPredict$FIRREratio == sexPredict$XISTratio
sexPredict$concordantLabel<-xPeaksWide$sex == sexPredict$FIRREratio
sexPredict$concordantLabel[sexPredict$concordantPredictions != TRUE]<-FALSE

xPeaksWide$qcpass2 <- sexPredict[match(xPeaksWide$sampleID, sexPredict$sampleID),]$concordantLabel
xPeaksWide <- xPeaksWide[complete.cases(xPeaksWide[ , c('qcpass2')]), ]

## Supplementary Figure 4.a: Results sex check ##
## x counts nearby peaks in FIRRE/XIST
## y total counts in peaks in Y chromosome 
## color by phenotype sex
## shape: predicted vs real sex

a1<-ggplot(xPeaksWide, aes(x=FIRRE, y=YTot, color=sex, shape = qcpass2)) + 
  geom_point() + 
  labs(
    x = "Counts in FIRRE associated peak",
    y = "Counts across all Y chr peaks", color = "Sex", title = "A", shape = "PASS QC"
  ) +scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16))+
   scale_y_continuous(
    labels = label_number(suffix = "k", scale = 1e-3))+
 pub_theme + geom_abline(intercept = 0, slope = thresF) 
  
a2 <- ggplot(xPeaksWide, aes(x=XIST, y=YTot, color=sex, shape = qcpass2)) + 
  geom_point() + 
  labs(
    x = "Counts in XIST associated peak",
    y = "Counts across all X chr peaks",color = "Sex", shape = "PASS QC"
  ) +scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16))+
   scale_y_continuous(labels = label_number(suffix = "k", scale = 1e-3))+
   pub_theme + geom_abline(intercept = 0, slope = thresX)

## Supplementary Figure 4.b: Results genotype check ##
## x number of reads
## y FREEMIX/CHIPMIX scores 
## color by pass/fail

genoFolder<-file.path(dir, "3_aligned/genotypeConcordance")

setwd(genoFolder)
outFiles<-list.files(genoFolder, pattern = "selfSM$")

genoScores<-NULL
for(each in outFiles){
	genoScores<-rbind(genoScores, cbind(sub("\\.selfSM", "", each), read.table(each)))
}

colnames(genoScores)<-c("sampleID", "SEQ_ID", "RG", "CHIP_ID", "nSNPS", "nREADS", "AVG_DP", "FREEMIX", "FREELK1", "FREELK0", "FREE_RH", "FREE_RA", "CHIPMIX","CHIPLK1","CHIPLK0","CHIP_RH","CHIP_RA","DPREF","RDPHET","RDPALT")

genoScores<-genoScores[match(pheno$sampleID, genoScores$sampleID),]
fraction <- pheno[match(genoScores$sampleID, pheno$sampleID),]$fraction
genoScores$fraction <- fraction
genoScores$cohort <- pheno[match(genoScores$sampleID, pheno$sampleID),]$cohort
genoScores$batch <- pheno[match(genoScores$sampleID, pheno$sampleID),]$sequencingBatch
genoScores$sampleCode <- pheno[match(genoScores$sampleID, pheno$sampleID),]$sampleCode
genoScores <- genoScores[complete.cases(genoScores[ , c('FREEMIX', 'CHIPMIX')]), ]
#genoScores$qcpass2 <- pheno[match(genoScores$sampleID, pheno$sampleID),]$qcpass2
#genoScores$qcpass2 <- as.factor(genoScores$qcpass2)

genoScores$qcpass <- ifelse(genoScores$FREEMIX < 0.1, "TRUE", "FALSE")
genoScores.plots.p1 <- genoScores                   
genoScores.plots.p1$qcpass <- genoScores.plots.p1$qcpass  

p1 <- ggplot(genoScores.plots.p1, aes(x=nREADS, y=FREEMIX, color=qcpass)) +
  geom_point(size=2) +
  geom_hline(yintercept=0.1) +
  scale_color_manual(values = col_vals_pass) +
  labs(x="Number of reads", y="FREEMIX", color="PASS QC", title = "B") +  
  scale_x_continuous(labels = label_number(suffix = "M", scale = 1e-5)) +
  pub_theme

genoScores$qcpass <- ifelse(genoScores$CHIPMIX < 0.1, "TRUE", "FALSE")
genoScores.plots.p2 <- genoScores
genoScores.plots.p2$qcpass <- genoScores.plots.p2$qcpass  

p2 <- ggplot(genoScores.plots.p2, aes(x=nREADS, y=CHIPMIX, color=qcpass)) +
  geom_point(size=2) +
  geom_hline(yintercept=0.1) +
  scale_color_manual(values = col_vals_pass) +
  labs(x="Number of reads", y="CHIPMIX", color="PASS QC") +
  scale_x_continuous(labels = label_number(suffix = "M", scale = 1e-5)) +
  pub_theme

## Combined figures ##
combined_fig5 <- (a1 + a2 + 
  plot_layout( guides = "collect", nrow = 1))& theme(
  legend.position = "bottom"
)

combined_fig52 <- (p1 + p2 + 
  plot_layout( guides = "collect", nrow = 1))& theme(
  legend.position = "bottom"
)
combFig5 <- combined_fig5 / combined_fig52 
 
ggsave(
  filename = paste0("SupplementaryFig_4.png"),
  plot = combFig5,
  width = 10,
  height = 10
)

ggsave(
  filename = paste0("SupplementaryFig_4.pdf"),
  plot = combined_fig52,
  width = 10,
  height = 10
)
