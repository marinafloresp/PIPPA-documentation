######################################
##       Supplementary Fig 2,3      ##
######################################

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


dfa <- data.frame(
  Samples = qcStats$sampleID,
  qcpass = qcStats$qcpass,
  reads = qcStats$NumReadsbeforeQC
) 

plot.a <- ggplot(dfa, aes( x = reads)) +
  geom_histogram(bins = 100,fill ='gray')+
  labs( x = expression("Number of reads (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", title = "A") +
  pub_theme +
  scale_x_continuous(
    labels = label_number(suffix = "", scale = 1e-6)
  )

dfb <- data.frame(
  Samples = qcStats$sampleID,
  reads = qcStats$UniqueNumReadsQC
)

plot.b <- ggplot(dfb,aes(x =  reads))+
geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% unique reads', y = "Number of samples", title = "B")+
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

qcStats$qcpass <- ifelse(qcStats$TotalNumReads > numReads*10^6, "TRUE", "FALSE")
plot.d <- ggplot(qcStats, aes(x = TotalNumReads,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=numReads*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x =  expression("Number of reads aligned (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", fill = "PASS",title = "D")+
  pub_theme +
	scale_x_continuous(labels = label_number(suffix = " ", scale = 1e-6))

plot.e <- ggplot(qcStats, aes(x =  PercentMTReads))+
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
plot.k <- plot.k + theme(legend.position = "bottom")


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
)
  )

ggsave(
  filename = paste0("SupplementaryFig_2_D1.png"),
  plot = combined_figall,
  width = 13,
  height = 17)
  
ggsave(
  filename = paste0("SupplementaryFig_2_D1.pdf"),
  plot = combined_figall,
  width = 13,
  height = 17
)



