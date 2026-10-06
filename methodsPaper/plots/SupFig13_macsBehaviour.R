#####################################
##       Supplementary Fig 13      ##
#####################################

library(kneedle)
library(segmented)
library(tidyr)
library(dplyr)

## Universal theme for plotting
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.key.width = unit(2.2, "cm")
  )

col_vals <- c("SOX10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")

######################
## gather peak sets ##
######################

groupsPeaks <-"/methodsPaper/data/relaxedPeakSets"
peakFiles2 <- c("NEUN_5e-2.narrowPeak.filt","SOX_5e-2.narrowPeak.filt",
                "neun_q5e-2.narrowPeak.filt","sox10_q5e-2.narrowPeak.filt"
)

allPeakFiles <- paste0(groupsPeaks,"/",peakFiles2)
peaks<-c()
for(file in allPeakFiles){
  peaks[[file]]<-makeGRangesFromDataFrame(read.table(file.path(file)), keep.extra.columns=TRUE, seqnames.field="V1",
                                          start.field="V2",
                                          end.field="V3",
                                          ignore.strand=TRUE)
}

names(peaks) <- gsub(groupsPeaks, "", names(peaks))
names(peaks) <- gsub(paste0(groupsPeaks, "/"), "", names(peaks))
names(peaks) <- c("NeuN+_Dataset1","SOX10+_Dataset1","NeuN+_Dataset2","SOX10+_Dataset2")

cf <- "SOX"
peakSets <- peaks$`SOX10+_Dataset1`
load(paste0(cf,"_5e-2_passAll.narrowPeak.rdata"))
df <-  as.data.frame(peak_set_group)
peakSets$PeakID2 <- paste0(peakSets$chr, "-",peakSets$start, "-", peakSets$end)
peakSets$n_samples <- df[match(peakSets$PeakID2,rownames(df)),]$n_samples
peaks.df <- peakSets
peakSets$fraction <- cf
df.list <- peaks.df

# ########################## Fig Sup 13.a ###########################################
# Violin plots of distribution of number of samples overlap with (y-axis)
# q-value threshold (x-axis). 
# coloured by cell-fraction (facets)

exponents <- seq(-10, -100, by = -5)

thresholds <- sort(c(0.05, 10^exponents))

breaks <- c(-Inf, thresholds)

labels <- as.character(thresholds)
#peaks.df <- peakSets
peaks.df$qValuegroup <- cut(
  peaks.df$qvals,
  breaks = breaks,
  labels = labels,
  include.lowest = TRUE
)
#peaks.df$qValuegroup[is.na(peaks.df$qValuegroup)] <- "1e-100"
peaks.df$threshold_log <- as.character(peaks.df$qValuegroup)
peaks.df$threshold_exp <- ifelse(
  peaks.df$threshold_log == "0.05",
  1.3,
  as.numeric(sub("1e-", "", peaks.df$threshold_log))
)
peaks.df$threshold_exp <- factor(
  peaks.df$threshold_exp,
  levels = c(1.3, seq(10, 100, 5))
)

p13a <- ggplot(peaks.df, aes(x = threshold_exp, y = n_samples))+ 
  geom_violin(scale = "width", fill = "#6B8E23")+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+
  scale_x_discrete(
    labels = c("1.3", seq(10, 100, 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+
  labs(y = "Number of samples overlap", title = "A") +
  pub_theme +
  theme(
    legend.position = "none")

# ########################## Fig Sup 13.b ###########################################
# Line graph of number of peaks (y-axis)
# q-value threshold (x-axis) 
# with one line per “peak size group”.  
# coloured by cell-fraction (facets)

exponents <- seq(-4, -100, by = -2)
thresholds <- sort(c(0.05,10^exponents))
labels <- as.character(thresholds)

breaks <- c(thresholds, -Inf)


peaks.df$qValuegroup <- cut(
  peaks.df$qvals,
  breaks = breaks,
  labels = labels,
  include.lowest = TRUE
)

peaks.df$qValuegroup <- factor(
  peaks.df$qValuegroup,
  levels = rev(labels)
)

peaks.df$group <- cut(
  peaks.df$width,
  breaks = c(0,400,700,1200,2500,4000,Inf),
  labels = c(
    "very narrow\n(0-400)",
    "narrow\n(400-700)",
    "medium\n(700-1200)",
    "broad\n(1200-2500)",
    "very broad\n(2500-4000)",
    "extreme\n(>4000)"
  ),
  include.lowest = TRUE
)

plot_df <- peaks.df %>% dplyr::count(qValuegroup, group)
plot_df$threshold_log <- as.character(plot_df$qValuegroup)
plot_df$threshold_log <- as.numeric(plot_df$threshold_log)
plot_df$threshold_log <- -log10(plot_df$threshold_log)

plot_df$threshold_exp <- ifelse(
  plot_df$threshold == "0.05",
  1.3,
  as.numeric(sub("1e-", "", plot_df$threshold_log))
)

peaks.df$group <- cut(
  peaks.df$width,
  breaks = c(0,400,700,1200,2500,4000,Inf),
  labels = c(
    "very_narrow\n(0-400)",
    "narrow\n(400-700)",
    "medium\n(700-1200)",
    "broad\n(1200-2500)",
    "very_broad\n(2500-4000)",
    "extreme\n(>4000)"
  ),
  include.lowest = TRUE
)
colnames(peaks.df)[13] <- "qvals2"
plot_df <- peaks.df %>% dplyr::count(qValuegroup, group)
plot_df$threshold_log <- as.character(plot_df$qValuegroup)
plot_df$threshold_exp <- ifelse(
  plot_df$threshold == "0.05",
  1.3,
  as.numeric(sub("1e-", "", plot_df$threshold_log))
)
p13b <-ggplot(plot_df, aes(x = threshold_log, y = n, group=group, color = group)) +
  geom_line(linewidth = 1)+
  labs(
    x = "Q value",
    y = "Number of peaks", title = "B", color = "Peak width\ngroup\n(bp)"
  ) +scale_color_manual(values = colorBlindGrey8)+
  scale_x_continuous(
    breaks = c(1.3, seq(10, 100, by = 5)),
    labels = c("1.3",seq(10, 100, by = 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+
  theme_minimal(base_size = 14) +
  pub_theme +
  theme(
    legend.position = "bottom")

# ########################## Fig Sup 13.c ###########################################
# Line graph of number of peaks (y-axis) 
# q-value threshold (x-axis) with
# one line per overlap with 0,1,2,3 samples etc.  
# plot 0 up to 10 individually, then 10-20, 20-30, or more than 30

peaks.df$n_samples_group <- cut(
  as.numeric(peaks.df$n_samples),
  breaks = c(0,1,2,3,6,10, 20, 30, Inf),
  labels = c("0","1","2","3-6","7-10","10–20", "20–30", ">30"),
  right = TRUE,
  include.lowest = TRUE
)

table(peaks.df$qValuegroup)   # q-value bins (factor)
table(peaks.df$n_samples_group)     # number of samples per peak
peaks.df$qValuegroup <- factor(peaks.df$qValuegroup)
peaks.df$n_samples_group <- factor(peaks.df$n_samples_group)
tab <- table(peaks.df$qValuegroup, peaks.df$n_samples_group)#, peaks.df$fraction)
df_plot <- as.data.frame(tab)
colnames(df_plot) <- c("qValuegroup", "n_samples",  "count")#"fraction",
df_plot$threshold_log <- as.character(df_plot$qValuegroup)
df_plot$threshold_exp <- ifelse(
  df_plot$threshold_log == "0.05",
  1.3,
  as.numeric(sub("1e-", "", df_plot$threshold_log))
)

p13c <-ggplot(df_plot, aes(x = threshold_exp, y = count, group=n_samples, color = n_samples)) +
  #geom_point() +
  geom_line(linewidth = 1)+
  scale_x_continuous(
    breaks = c(-log10(0.05), seq(10, 100, by = 5)),
    labels = c("1.3", seq(10, 100, by = 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+scale_color_manual(values = colorBlindGrey8)+ 
  labs(
    y = "Number of peaks", title = "C", color = "Number of samples\noverlap"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

design <- "AAAAAA
           AAAAAA
           BBBCCC
           BBBCCC"

supFig13 <- free(p13a) + p13b + p13c +  plot_layout(design = design)

ggsave(
  filename = paste0("SupplementaryFig_13_D1_SOX.png"),
  plot = supFig13,
  width = 17,
  height = 12
)
ggsave(
  filename = paste0("SupplementaryFig_13_D1_SOX.pdf"),
  plot = supFig13,
  width = 17,
  height = 12
)
