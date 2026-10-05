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
qcStats <- read.csv("qcStage1Results_D2.csv") 

# ########################## Sup Fig 3.a #################
# Histogram of number of reads before QC per sample 

dfa <- data.frame(
  Samples = qcStats$Sample,
  reads = qcStats$Total.Sequences
) 

plot.a <- ggplot(dfa, aes( x = reads)) +
  geom_histogram(bins = 100,fill ='gray')+
  labs( x = expression("Number of reads (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", title = "A") +
  pub_theme +
  scale_x_continuous(
    labels = label_number(suffix = " ", scale = 1e-6)
  )

# ########################## Sup Fig 3.b #################
# Histogram of number of reads before QC per sample 

align_f <- as.data.frame(qcStats)
align_f$qcpass<- ifelse(align_f$overall_alignment_rate > 80, "TRUE", "FALSE")

plot.b <- ggplot(align_f, aes(x = overall_alignment_rate,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=80,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +
  labs(x  = '% Alignment rate', y = "Number of samples", title = "B", 
       fill = "PASS")+
  pub_theme +
  theme(legend.position = "none")
numReads <- 50
qcStats$qcpass <- ifelse(qcStats$aligned_total > numReads*10^6, "TRUE", "FALSE")
plot.c <- ggplot(qcStats, aes(x = aligned_total,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=numReads*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x = expression("Number of reads aligned (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", 
       fill = "PASS",title = "C")+
  pub_theme +
  scale_x_continuous(labels = label_number(suffix = "", scale = 1e-6))

plot.d <- ggplot(qcStats, aes(x =  PercentMTReads))+
  geom_histogram(bins = 100,fill ='gray')+
  labs(x  = '% Mitochondrial reads', y = "Number of samples", title = "D")+
  pub_theme 


eMetrics_f <- qcStats
eMetrics_f$qcpass <- ifelse(qcStats$PBC1 > 0.7, "TRUE", "FALSE")

eMetrics_g <- qcStats
eMetrics_g$qcpass <- ifelse(qcStats$PBC2 > 1, "TRUE", "FALSE")

eMetrics_h <- qcStats
eMetrics_h$qcpass <- ifelse(qcStats$NRF > 0.7, "TRUE", "FALSE")

plot.d <- ggplot(eMetrics_h, aes(x = NRF, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NRF', y = "Number of samples", title = "D", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.7, color = "red")+ theme(
    legend.position = "none"
  )

plot.e<- ggplot(eMetrics_f, aes(x = PBC1, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC1', y = "Number of samples", title = "E", fill = "PASS Stage 1 QC") +
  pub_theme + geom_vline(xintercept = 0.7, color = "red")

plot.f <- ggplot(eMetrics_g, aes(x = PBC2, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC2', y = "Number of samples", title = "F", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 1, color = "red")


propNuc_j <-qcStats
propNuc_j$qcpass <- ifelse(qcStats$NFR > 0.25, "TRUE", "FALSE")

propNuc_k <-qcStats
propNuc_k$qcpass <- ifelse(qcStats$Mono > 0.2, "TRUE", "FALSE")

plot.g <- ggplot(propNuc_j, aes(x = NFR, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NFR fragments', y = "Number of samples", title = "G", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.25, color = "red") + 
  theme(legend.position = "none")

plot.h <- ggplot(propNuc_k, aes(x = Mono, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'MNC fragments', y = "Number of samples", title = "H", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.2, color = "red")

plot.c <- plot.c + theme(legend.position = "none")
plot.f <- plot.f + theme(legend.position = "none")
plot.e <- plot.e + theme(legend.position = "none")

plot.d <- plot.d + theme(legend.position = "none")
plot.g <- plot.g + theme(legend.position = "none")
plot.h <- plot.h + theme(legend.position = "none")

legend_shared <- cowplot::get_legend(
  plot.e +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
)
design = "
AAABBB
AAABBB
CCCDDD
CCCDDD
EEEFFF
EEEFFF
GGGHHH
GGGHHH
..LL..
"
combined_figall <- (plot.a + plot.b + plot.c +
                      plot.d + plot.e + plot.f +
                      plot.g + plot.h +  legend_shared + 
                      plot_layout(ncol = 2, heights = c(1,1,1,1,0.25))+ 
plot_layout(guides = "collect"             # <-- collect all legends into guide_area()
) & theme(legend.position = "none"))
