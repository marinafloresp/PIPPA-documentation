#########################
##       Figure 3      ##
#########################

library(tidyr)
library(dplyr)
library(GenomicRanges)

groupsPeaks <- "/dir/to/peak/sets/"

format_q_label <- function(q) {
  
  exponent <- round(log10(q))
  
  mantissa <- q / (10^exponent)
  
  nearest <- c(1, 5)[
    which.min(abs(c(1, 5) - mantissa))
  ]
  
  paste0(nearest, "e", exponent)
}

## Universal theme for plotting
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 24),
    legend.key.width = unit(2.2, "cm")
  )
#color blind friendly palette for multiple colours
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")


col_vals <- c("SOX10+" = "#6B8E23","NeuN+" = "hotpink")

######################
## Gather peak sets ##
######################
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
lapply(peaks, nrow)
width <- lapply(peaks, width)
totalwidth <- lapply(width, sum)

names(peaks) <- c("NeuN+_Dataset1","SOX10+_Dataset1","NeuN+_Dataset2","SOX10+_Dataset2")

peaks <- lapply(peaks, function(df) {
  df <- as.data.frame(df)
  colnames(df) <-c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal", "pvalue", "qvalue", "score2")
  df
})

peaks <- lapply(peaks, function(df) {
  df$qvals <- 10^(-df$qvalue)
  df
})

exponents <- seq(-2, -100, by = -2)
thresholds <- c(0.05,10^exponents)
order_levels.og <- as.character(thresholds)


counts_list <- lapply(peaks, function(gr) {
  qvals <- gr$qvals
  sapply(thresholds, function(t) sum(qvals < t, na.rm = TRUE))
})


df <- do.call(rbind, lapply(seq_along(peaks), function(i) {
  data.frame(
    set = names(peaks)[i],
    threshold = thresholds,
    count = counts_list[[i]]
  )
}))

df[c("Cell-fraction", "Dataset")] <- do.call(
  rbind,
  strsplit(df$set, "_")
)
df$threshold_log <- -log10(df$threshold)
p3a <- ggplot(df, aes(
  x = threshold_log,
  y = count,
  color = `Cell-fraction`,
  linetype = Dataset,
  group = interaction(`Cell-fraction`, Dataset)
)) + scale_color_manual(values= col_vals)+
  scale_x_continuous(
    breaks = c(1.3, seq(10, 100, by = 5)),
    labels = c("1.3", seq(10, 100, by = 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+
  geom_line(linewidth = 1) +
  labs(
    x = "Q value",
    y = "Number of peaks",
    title = "A", color = "Cell type"
  ) +
  # scale_x_discrete(limits = order_levels.og) +
  theme_minimal(base_size = 14) +
  pub_theme +
  theme(
    legend.position = "bottom")

# ########################## Fig 3.b ###########################################
# Line graph of number of peaks (y-axis)
# q-value threshold (x-axis) 
# with one line per “peak size group”.  
peaks.df <- peaks$`NeuN+_Dataset1`
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
peaks.df %>% dplyr::count(qValuegroup)

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

p3b <-ggplot(plot_df, aes(x = threshold_log, y = n, group=group, color = group)) +
  #  geom_point() +
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

# ########################## Fig 3.c ###########################################
# Line graph of number of peaks (y-axis) 
# q-value threshold (x-axis) with
# one line per overlap with 0,1,2,3 samples etc.  
# plot 0 up to 10 individually, then 10-20, 20-30, or more than 30

load("NEUN_5e-2_passAll.narrowPeak.rdata")
peak_set_group.df <- as.data.frame(peak_set_group)
peaks.df$PeakID2 <- paste0(peaks.df$chr, "-",peaks.df$start, "-", peaks.df$end)
peaks.df$n_samples <- peak_set_group.df[match(peaks.df$PeakID2,rownames(peak_set_group.df)),]$n_samples

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
tab <- table(peaks.df$qValuegroup, peaks.df$n_samples_group)
df_plot <- as.data.frame(tab)
colnames(df_plot) <- c("qValuegroup", "n_samples", "count")
df_plot$threshold_log <- as.character(df_plot$qValuegroup)
df_plot$threshold_log <- as.numeric(df_plot$threshold_log)
df_plot$threshold_log <- -log10(df_plot$threshold_log)

df_plot$threshold_exp <- ifelse(
  df_plot$threshold_log == "0.05",
  1.3,
  as.numeric(sub("1e-", "", df_plot$threshold_log))
)

p3c <-ggplot(df_plot, aes(x = threshold_log, y = count, group=n_samples, color = n_samples)) +

  geom_line(linewidth = 1)+
  scale_x_continuous(
    breaks = c(1.3, seq(10, 100, by = 5)),
    labels = c("1.3", seq(10, 100, by = 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+scale_color_manual(values = colorBlindGrey8)+
  labs(
    y = "Number of peaks", title = "C", color = "Number of samples\noverlap"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom",
  )

# ########################## Fig 3d ###########################################
# Violin plots of distribution of number of samples overlap with (y-axis)
# q-value threshold (x-axis). 


exponents <- seq(-10, -100, by = -5)
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

p3d<- ggplot(peaks.df, aes(x = threshold_exp, y = n_samples))+ 
  geom_violin(scale = "width", fill ="hotpink")+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+
  scale_x_discrete(
    labels = c("1.3", seq(10, 100, 5)),
    name = expression("Q value threshold (" * 10^{-n} * ")")
  )+
  labs( y = "Samples overlap", title = "D") +
  pub_theme +
  theme(
    legend.position = "none")

## combine plots for figure 3 ##
design <- "AAABBB
           AAABBB
           CCCDDD
           CCCDDD"

figure3 <- p3a + p3b + p3c+ p3d +  plot_layout(design = design)


## Save figure ##
ggsave(
  filename = paste0("Figure_3.pdf"),
  plot = figure3,
  width = 22,
  height = 12
)

ggsave(
  filename = paste0("Figure_3.png"),
  plot = figure3,
  width = 22,
  height = 12
)
