###################################
##       Supplementary Fig 3     ##
###################################

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

p3a <- ggplot(dfa, aes( x = reads)) +
  geom_histogram(bins = 100,fill ='gray')+
  labs( x = expression("Number of reads (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", title = "A") +
  pub_theme +
  scale_x_continuous(
    labels = label_number(suffix = " ", scale = 1e-6)
  )

# ########################## Sup Fig 3.b #################
# Histogram of Alignment rate metric per sample
# Samples coloured by pass/fail based on threshold (>80%), shown as vertical red line

align_f <- as.data.frame(qcStats)
align_f$qcpass<- ifelse(align_f$overall_alignment_rate > 80, "TRUE", "FALSE")

p3b <- ggplot(align_f, aes(x = overall_alignment_rate,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=80,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +
  labs(x  = '% Alignment rate', y = "Number of samples", title = "B", 
       fill = "PASS")+
  pub_theme +
  theme(legend.position = "none")

# ########################## Sup Fig 3.c #################
# Histogram of number of reads aligned (after QC) per sample
# Samples coloured by pass/fail based on threshold (>50M reads), shown as vertical red line

qcStats$qcpass <- ifelse(qcStats$aligned_total > numReads*10^6, "TRUE", "FALSE")
p3c <- ggplot(qcStats, aes(x = aligned_total,fill = qcpass))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept=numReads*10^6,
             color ='red')+
  scale_fill_manual(values = col_vals_pass) +           
  labs(x = expression("Number of reads aligned (" ~ x ~ 10^6 ~ ")"), y = "Number of samples", 
       fill = "PASS",title = "C")+
  pub_theme +
  scale_x_continuous(labels = label_number(suffix = "", scale = 1e-6))

# ########################## Sup Fig 3.d #################
# Histogram of NRF metric per sample
# Samples coloured by pass/fail based on threshold (>0.7), shown as vertical red line

eMetrics_h <- qcStats
eMetrics_h$qcpass <- ifelse(eMetrics_h$NRF > nrfThres, "TRUE", "FALSE")

p3d <- ggplot(eMetrics_h, aes(x = NRF, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NRF', y = "Number of samples", title = "D", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.7, color = "red")+ theme(
    legend.position = "none"
  )

# ########################## Sup Fig 3.e #################
# Histogram of PBC1 metric per sample
# Samples coloured by pass/fail based on threshold (>0.7), shown as vertical red line

eMetrics_f <- qcStats
eMetrics_f$qcpass <- ifelse(eMetrics_f$PBC1 > pbc1Thres, "TRUE", "FALSE")

p3e<- ggplot(eMetrics_f, aes(x = PBC1, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC1', y = "Number of samples", title = "E", fill = "PASS Stage 1 QC") +
  pub_theme + geom_vline(xintercept = 0.7, color = "red")

# ########################## Sup Fig 3.f #################
# Histogram of PBC2 metric per sample
# Samples coloured by pass/fail based on threshold (>1), shown as vertical red line

eMetrics_g <- qcStats
eMetrics_g$qcpass <- ifelse(eMetrics_g$PBC2 > pbc2Thres, "TRUE", "FALSE")

p3f <- ggplot(eMetrics_g, aes(x = PBC2, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'PBC2', y = "Number of samples", title = "F", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 1, color = "red")

# ########################## Sup Fig 3.g #################
# Histogram of %NFR metric per sample
# Samples coloured by pass/fail based on threshold (>0.25), shown as vertical red line

propNuc_j <-qcStats
propNuc_j$qcpass <- ifelse(qcStats$NFR > nucfThres, "TRUE", "FALSE")

p3g <- ggplot(propNuc_j, aes(x = NFR, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'NFR fragments', y = "Number of samples", title = "G", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.25, color = "red") + 
  theme(legend.position = "none")

# ########################## Sup Fig 3.h #################
# Histogram of %MNC metric per sample
# Samples coloured by pass/fail based on threshold (>0.2), shown as vertical red line

propNuc_k <-qcStats
propNuc_k$qcpass <- ifelse(qcStats$monoThres > 0.2, "TRUE", "FALSE")

p3h <- ggplot(propNuc_k, aes(x = Mono, fill = qcpass)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col_vals_pass) +
  labs(x = 'MNC fragments', y = "Number of samples", title = "H", fill = "PASS") +
  pub_theme + geom_vline(xintercept = 0.2, color = "red")

## Combine all plots for figure ##
p3c <- p3c + theme(legend.position = "none")
p3f <- p3f + theme(legend.position = "none")
p3e <- p3e + theme(legend.position = "none")
p3d <- p3d + theme(legend.position = "none")
p3g <- p3g + theme(legend.position = "none")
p3h <- p3h + theme(legend.position = "none")
legend_shared <- cowplot::get_legend(
  p3e +
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
supFig3 <- (p3a + p3b + p3c +
                      p3d + p3e + p3f +
                      p3g + p3h +  legend_shared + 
                      plot_layout(ncol = 2, heights = c(1,1,1,1,0.25))+ 
plot_layout(guides = "collect"             
) & theme(legend.position = "none"))

ggsave(
  filename = paste0("SupplementaryFig_3_D2.png"),
  plot = supFig3,
  width = 13,
  height = 17)
  
ggsave(
  filename = paste0("SupplementaryFig_3_D2.pdf"),
  plot = supFig3,
  width = 13,
  height = 17
)
