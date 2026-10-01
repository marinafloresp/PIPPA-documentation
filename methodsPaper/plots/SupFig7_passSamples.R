#######################################
##       Supplementary Figure 7      ##
#######################################

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
  library(plotly)
  library(tidyr)
  library(pheatmap)
  library(cowplot)
  library(readr)
  library(patchwork)
}))

## settings for plotting ##
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")
pass_fail_cols <- c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
)

#######################
## Gather QC results ##
#######################

qcStatsD2 <- read_csv("./methodsPaper/data/QC1_D2.csv")
qcStats2D2 <- read_csv("./methodsPaper/data/MS_alignmentRate.csv")
qcStatsD2 <- qcStatsD2[,-1]
qcStats2D2 <- qcStats2D2[,-1]
qc_cols <- c("PASS_AlignmentRate", "PASS_AlignedTotal", "PASS_NRF", "PASS_PBC1","PASS_PBC2","PASS_NFR","PASS_Mono" )
qcStatsD2$PASSALL <-  as.factor(rowSums(qcStatsD2[, qc_cols]) == length(qc_cols))
qcStatsD2$Ratio <- qcStatsD2$NFR/qcStatsD2$Mono
qcStatsD2$aligned_total_reads <- qcStatsD2$aligned_total*2
qcStatsD2[order(qcStatsD2$aligned_total, decreasing = TRUE),]$Sample

passAllQCStatus <- as.data.frame(read_csv("./methodsPaper/passAllQCStatus.csv"))
qc1StatsD1 <- read.table("./methodsPaper/qc1Stats.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
qc1StatsD1$PASSALL <- passAllQCStatus[match(qc1StatsD1$sampleID,passAllQCStatus$sampleID),]$PASSALL
qc1StatsD1$Ratio <- qc1StatsD1$NFR/qc1StatsD1$MNC

## plot number of reads
D1.sub <- qc1StatsD1[,c("sampleID","TotalNumReads","PASSALL")]
D1.sub$dataset <- "Dataset 1"
D2.sub <- qcStatsD2[,c("Sample","aligned_total","PASSALL")]
D2.sub$dataset <- "Dataset 2"
colnames(D2.sub) <- colnames(D1.sub)
reads.df <- rbind(D1.sub,D2.sub)

reads.df$TotalNumReads <- reads.df$TotalNumReads/1000000

# ########################## Fig Sup 7.a ##################
# Violin plot of number of reads (y-axis)
# by dataset (x-axis) 
# coloured by whether sample passed/failed QC

p7a <- ggplot(
  reads.df,
  aes(
    x = dataset,
    y = TotalNumReads, fill = PASSALL)
) +
  geom_violin(scale = "width",
              position = position_dodge(width = 0.9),
              trim = FALSE) +
  geom_boxplot(width = 0.15,
               position = position_dodge(width = 0.9),
               outlier.shape = NA)+
  stat_summary(
    fun = median,position = position_dodge(width = 0.9),
    geom = "point",
    size = 2.5,
    color = "white"
  )+
  scale_fill_manual(values = pass_fail_cols) +
  labs(
    x = "",
    y =expression("Number of reads (" ~ x ~ 10^6 ~ ")"), title = "A"
  ) +
  pub_theme +  theme_bw(12)+
  theme(axis.title = element_text(size = 14), legend.position = "none",
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 24,face = "bold"),
        plot.subtitle = element_text(size = 22,
                                     face = "bold",
                                     hjust = 0.5))

# ########################## Fig Sup 7.b ##################
# Violin plot of NFR to MNC ratio (y-axis)
# by dataset (x-axis) 
# coloured by whether sample passed/failed QC

D1.sub <- qc1StatsD1[,c("sampleID","Ratio","PASSALL")]
D1.sub$dataset <- "Dataset 1"
D2.sub <- qcStatsD2[,c("Sample","Ratio","PASSALL")]
D2.sub$dataset <- "Dataset 2"
colnames(D2.sub) <- colnames(D1.sub)
reads.df <- rbind(D1.sub,D2.sub)

p7b <- ggplot(
  reads.df,
  aes(
    x = dataset,
    y = Ratio, fill = PASSALL)
) +
  geom_violin(scale = "width",
              position = position_dodge(width = 0.9),
              trim = FALSE) +
  geom_boxplot(width = 0.15,
               position = position_dodge(width = 0.9),
               outlier.shape = NA)+
  stat_summary(
    fun = median,position = position_dodge(width = 0.9),
    geom = "point",
    size = 2.5,
    color = "white"
  )+
  scale_fill_manual(values = pass_fail_cols) +
  labs(
    x = "",
    y = "NFR/MNC", title = "B", fill = "Pass QC"
  ) +
  pub_theme +  theme_bw(12)+
  theme(axis.title = element_text(size = 14), legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 24,face = "bold"),
        plot.subtitle = element_text(size = 22,
                                     face = "bold",
                                     hjust = 0.5))
p7b

# ########################## Fig Sup 7.c ##################
# Violin plot of TSS enrichment scores (y-axis)
# by dataset (x-axis) 
# coloured by whether sample passed/failed QC

D1.sub <- qc1StatsD1[,c("sampleID","TSS_score","PASSALL")]
D1.sub$dataset <- "Dataset 1"
p7c <- ggplot(
  D1.sub,
  aes(
    x = dataset,
    y = TSS_score, fill = PASSALL)
) +
  geom_violin(scale = "width",
              position = position_dodge(width = 0.9),
              trim = FALSE) +
  geom_boxplot(width = 0.15,
               position = position_dodge(width = 0.9),
               outlier.shape = NA)+
  stat_summary(
    fun = median,position = position_dodge(width = 0.9),
    geom = "point",
    size = 2.5,
    color = "white"
  )+
  scale_fill_manual(values = pass_fail_cols) +
  labs(
    x = "",
    y = "TSS enrichment score", title = "C"
  ) +
  pub_theme +  theme_bw(12)+
  theme(axis.title = element_text(size = 14), legend.position = "none",
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 24,face = "bold"),
        plot.subtitle = element_text(size = 22,
                                     face = "bold",
                                     hjust = 0.5))
p7c

## combine plots for figure ##
figSup7 <- p7a | p7b | p7c + plot_layout(
  nrow = 1,
  guides = "collect"
) 

ggsave(
  filename = "SupplementaryFig_7.pdf",
  plot = figSup7,
  width = 15,
  height = 6
)

ggsave(
  filename = "SupplementaryFig_7.png",
  plot = figSup7,
  width = 15,
  height = 6
)

