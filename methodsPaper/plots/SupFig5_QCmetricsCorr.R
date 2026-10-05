#####################################
##       Supplementary Fig 5       ##
#####################################

suppressWarnings(suppressPackageStartupMessages({
  library(corrplot)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(cowplot)
  library(patchwork)
  library(scales)
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

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

## Gather all QC Stage 1 results ##
qcStats <- read.csv(file.path(qc1stats),header = TRUE, sep = ',', stringsAsFactors = FALSE)
qcStats$qcpass <- QCPASS.df[match(qcStats$sampleID, QCPASS.df$sampleID),]$QCPASS

# ########################## Sup Fig 5 ###########################################
# Correlation heatmap of QC metrics refined to those that are used to filter samples and QCPASS indicator variable.

corMergeStats<-cbind(qcStats$TotalNumReads,
                     qcStats$NRF,
                     qcStats$PBC1,
                     qcStats$PBC2,
                     qcStats$NFR,
                     qcStats$MNC,
                     qcStats$AlignmentRate)
head(corMergeStats)                
colnames(corMergeStats)<- c('Total\nreads', 'NRF', 'PBC1', 'PBC2', "%NFR", "%MNC", 'Alingment\nrate')

corMergeStats<- corMergeStats[, colSums(corMergeStats, na.rm = TRUE) != 0]

png("SupplementaryFig_5_D1.png", width = 2100, height = 2000, res = 300)

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

pdf("SupplementaryFig_5_D1.pdf",width = 11, height = 10)

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
