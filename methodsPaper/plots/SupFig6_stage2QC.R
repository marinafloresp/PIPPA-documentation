#####################################
##       Supplementary Fig 6       ##
#####################################
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

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

## settings for plotting ##
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
## colours for plotting
col_vals_pass <-c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
)

col_vals_sex <- c(
  "Female" = "#CC79A7",  # reddish purple
  "Male"   = "#0072B2"   # blue
)

## gather data ##
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
pheno$fraction[pheno$fraction == "SOX"] <- "SOX10+"
pheno$fraction[pheno$fraction == "NEUN"] <- "NeuN+"
pheno$fraction[pheno$fraction == "IRF8"] <- "IRF8+"
pheno$sequencingBatch <- as.factor(pheno$sequencingBatch)

qcPass1 <- read.csv(qc1pass, stringsAsFactors = FALSE, strip.white = TRUE, header = TRUE)
qcPass1$QCPASS <- as.factor(qcPass1$QCPASS)

qcPass2 <- read.csv(qc2pass, stringsAsFactors = FALSE, strip.white = TRUE, header = TRUE)
qcPass2$QCPASS <- as.factor(qcPass2$QCPASS)

pheno$qcpass1 <- qcPass1[match(pheno$sampleID,qcPass1$sampleID),]$QCPASS
pheno$qcpass2 <- qcPass2[match(pheno$sampleID,qcPass2$sampleID),]$QCPASS

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

xPeaksWide$sampleID<-gsub(".chrX.bam", "", xPeaksWide$sampleID, fixed = TRUE)
yPeaksWide$sampleID<-gsub(".chrY.bam", "", yPeaksWide$sampleID,fixed = TRUE)

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

xPeaksWide$sex[xPeaksWide$sex == "F"] <- "Female"
xPeaksWide$sex[xPeaksWide$sex == "M"] <- "Male"

# ########################## Supplementary Fig 6.a ###########################################
# Dot point for counts in peaks nearby genes XIST/FIRRE
# coloured by phenotype sex group
# shape is whether pass Stage 2 QC or not

a1<-ggplot(xPeaksWide, aes(x=FIRRE, y=YTot, color=sex, shape = qcpass2)) + 
  geom_point() + 
  labs(
    x = "Counts in FIRRE associated peak",
    y = "Counts across all Y chr peaks", color = "Sex", title = "A", shape = "PASS"
  ) +scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16))+
  scale_color_manual(values = col_vals_sex)+
   scale_y_continuous(
    labels = label_number(suffix = "k", scale = 1e-3))+
    guides(
    color = guide_legend(
      override.aes = list(size = 4)
    ),
    shape = guide_legend(
      override.aes = list(size = 4)
    )
  )+
 pub_theme + geom_abline(intercept = 0, slope = thresF,color = "red") 
  
a2 <- ggplot(xPeaksWide, aes(x=XIST, y=YTot, color=sex, shape = qcpass2)) + 
  geom_point() + 
  labs(
    x = "Counts in XIST associated peak",
    y = "Counts across all X chr peaks",color = "Sex", shape = "PASS"
  ) +scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16))+
  guides(
    color = guide_legend(
      override.aes = list(size = 4)
    ),
    shape = guide_legend(
      override.aes = list(size = 4)
    )
  )+
  scale_color_manual(values = col_vals_sex)+
   scale_y_continuous(labels = label_number(suffix = "k", scale = 1e-3))+
   pub_theme + geom_abline(intercept = 0, slope = thresX, color = "red")

# ########################## Supplementary Fig 6.b ###########################################
# Dot point for genotype check results
# FREEMIX value (x-axis)
# CHIPMIX value (y-axis)
# coloured by pass/fail genotype check

## Genotype check

genoFolder<-file.path(genoFolder)

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

genoScores$qcpassCHIPMIX <- genoScores$CHIPMIX < 0.1
genoScores$qcpassFREEMIX <- genoScores$FREEMIX < 0.1
genoScores$qcpass <- genoScores$qcpassFREEMIX & genoScores$qcpassCHIPMIX
genoScores.plots.p3 <- genoScores

p3 <- ggplot(genoScores.plots.p3, aes(x=FREEMIX, y=CHIPMIX, color=qcpass)) +
  geom_point(size=2) +
  geom_hline(yintercept=0.1,color = "red") +
  geom_vline(xintercept=0.1,color = "red") +
  scale_color_manual(values = col_vals_pass) +
  labs(x="FREEMIX", y="CHIPMIX", title = "B",color="PASS")  +
  pub_theme
 
combined_fig5 <- (a1 + a2 + 
  plot_layout( guides = "collect", nrow = 1))& theme(
  legend.position = "bottom"
)

combined_fig52 <- (p3 + 
  plot_layout( guides = "collect", nrow = 1))& theme(
  legend.position = "bottom"
)

design <- "AAAAAA
           AAAAAA
           #BBBB#
           #BBBB#
           "
           
combSupFig6 <- combined_fig5 / free(combined_fig52) + plot_layout( design = design, height= c(1,2))
 
ggsave(
  filename = paste0("SupplementaryFig_6.png"),
  plot = combSupFig6,
  width = 10,
  height = 12
)

ggsave(
  filename = paste0("SupplementaryFig_6.pdf"),
  plot = combSupFig6,
  width = 10,
  height = 12
)
