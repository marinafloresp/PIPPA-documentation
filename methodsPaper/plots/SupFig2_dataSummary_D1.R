####################################
##       Supplementary Fig 2      ##
####################################

suppressWarnings(suppressPackageStartupMessages({
  library(corrplot)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(cowplot)
  library(patchwork)
  library(scales)
}))

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

## Universal theme for plotting ## 

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

col_vals_pass <- c(
  "FALSE" = "#D55E00",  # vermillion (fail)
  "TRUE"  = "#009E73"   # bluish green (pass)
)

## Load file with QC Stage 1 results ##
qcStats <- read.csv("qcStage1Results.csv") 

# ########################## Sup Fig 2.a #################
# Histogram of number of reads before QC per sample 

dfa <- data.frame(
  Samples = qcStats$sampleID,
  qcpass = qcStats$qcpass,
  reads = qcStats$NumReadsbeforeQC
) 

p2a <- ggplot(dfa, aes( x = reads)) +
  geom_histogram(bins = 100,fill ='gray')+
  labs( x = expression("Number of reads (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", title = "A") +
  pub_theme +
  scale_x_continuous(
    labels = label_number(suffix = "", scale = 1e-6)
  )

# ########################## Sup Fig 2.b #################
# Histogram of % unique reads before QC per sample 

dfb <- data.frame(
  Samples = qcStats$sampleID,
  reads = qcStats$UniqueNumReadsQC
)

p2b <- ggplot(dfb,aes(x =  reads))+
geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% unique reads', y = "Number of samples", title = "B")+
  pub_theme 

# ########################## Sup Fig 2.c #################
# Histogram of Alignment rate metric per sample
# Samples coloured by pass/fail based on threshold (>80%), shown as vertical red line

align_f <- qcStats
align_f$qcpass<- ifelse(align_f$AlignmentRate > alignThres, "TRUE", "FALSE")

p2c <- ggplot(align_f, aes(x = AlignmentRate,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=80, 
             color ='red')+
             scale_fill_manual(values = col_vals_pass) +
  labs(x  = '% Alignment rate', y = "Number of samples", title = "C", fill = "PASS")+
  pub_theme 

# ########################## Sup Fig 2.d #################
# Histogram of number of reads aligned (after QC) per sample
# Samples coloured by pass/fail based on threshold (>50M reads), shown as vertical red line

qcStats$qcpass <- ifelse(qcStats$TotalNumReads > numReads*10^6, "TRUE", "FALSE")
p2d <- ggplot(qcStats, aes(x = TotalNumReads,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=numReads*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x =  expression("Number of reads aligned (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", fill = "PASS",title = "D")+
  pub_theme +
	scale_x_continuous(labels = label_number(suffix = " ", scale = 1e-6))

# ########################## Sup Fig 2.e #################
# Histogram of % of mitochondrial reads per sample

p2e <- ggplot(qcStats, aes(x =  PercentMTReads))+
  geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% MT reads', y = "Number of samples", title = "E")+
  pub_theme 

# ########################## Sup Fig 2.f #################
# Histogram of TSS enrichment score per sample

p2f<- ggplot(qcStats, aes(x  =TSS_score)) + 
  geom_histogram(bins = 100,fill ='gray')+
  labs(x =  'TSS enrichment score', y = "Number of samples", title = "F")+ 
  pub_theme 

# ########################## Sup Fig 2.g #################
# Histogram of PBC1 metric per sample
# Samples coloured by pass/fail based on threshold (>0.7), shown as vertical red line

eMetrics_f <- qcStats
eMetrics_f$qcpass <- ifelse(eMetrics_f$PBC1 > pbc1Thres, "TRUE", "FALSE")

p2g <- ggplot(eMetrics_f, aes(x = PBC1, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC1', y = "Number of samples", title = "G", fill = "PASS") +
  pub_theme + geom_vline(xintercept = pbc1Thres, color = "red")

# ########################## Sup Fig 2.h #################
# Histogram of PBC2 metric per sample
# Samples coloured by pass/fail based on threshold (>1), shown as vertical red line

eMetrics_g <- qcStats
eMetrics_g$qcpass <- ifelse(eMetrics_g$PBC2 > pbc2Thres, "TRUE", "FALSE")

p2h <- ggplot(eMetrics_g, aes(x = PBC2, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC2', y = "Number of samples", title = "H", fill = "PASS") +
  pub_theme + geom_vline(xintercept = pbc2Thres, color = "red")

# ########################## Sup Fig 2.i #################
# Histogram of NRF metric per sample
# Samples coloured by pass/fail based on threshold (>0.7), shown as vertical red line

eMetrics_h <- qcStats
eMetrics_h$qcpass <- ifelse(eMetrics_h$NRF > nrfThres, "TRUE", "FALSE")

p2i <- ggplot(eMetrics_h, aes(x = NRF, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NRF', y = "Number of samples", title = "I", fill = "PASS") +
  pub_theme + geom_vline(xintercept = nrfThres, color = "red")+ theme(
  legend.position = "none"
)

# ########################## Sup Fig 2.j #################
# Histogram of %NFR metric per sample
# Samples coloured by pass/fail based on threshold (>0.25), shown as vertical red line

propNuc_j <-qcStats
propNuc_j$qcpass <- ifelse(propNuc_j$NFR > nucfThres, "TRUE", "FALSE")

p2j <- ggplot(propNuc_j, aes(x = NFR, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NFR fragments', y = "Number of samples", title = "J", fill = "PASS") +
  pub_theme + geom_vline(xintercept = nucfThres, color = "red") + 
  theme(legend.position = "none")

# ########################## Sup Fig 2.k #################
# Histogram of %MNC metric per sample
# Samples coloured by pass/fail based on threshold (>0.2), shown as vertical red line

propNuc_k <-qcStats
propNuc_k$qcpass <- ifelse(propNuc_k$MNC > monoThres, "TRUE", "FALSE")

p2k <- ggplot(propNuc_k, aes(x = MNC, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'MNC fragments', y = "Number of samples", title = "K", fill = "PASS") +
  pub_theme + geom_vline(xintercept = monoThres, color = "red")

## combine all plots for figure ##

p2c <- p2c + theme(legend.position = "none")
p2d <- p2d + theme(legend.position = "none")
p2e<- p2e + theme(legend.position = "none")
p2f <- p2f + theme(legend.position = "none")
p2g <- p2g + theme(legend.position = "none")
p2h <- p2h + theme(legend.position = "none")
p2i <- p2i + theme(legend.position = "none")
p2j <- p2j + theme(legend.position = "none")
p2k <- p2k + theme(legend.position = "bottom")

supFig2 <- (p2a + p2b + p2c +
  p2d + p2e + p2f +
  p2h + p2i + p2f + 
  p2j + p2k + 
  plot_layout( design = "
AAABBB
AAABBB
CCCDDD
CCCDDD
EEEFFF
EEEFFF
GGHHII
GGHHII
JJJKKK
JJJKKK
"
)
  )

ggsave(
  filename = paste0("SupplementaryFig_2_D1.png"),
  plot = supFig2,
  width = 13,
  height = 17)
  
ggsave(
  filename = paste0("SupplementaryFig_2_D1.pdf"),
  plot = supFig2,
  width = 13,
  height = 17
)
