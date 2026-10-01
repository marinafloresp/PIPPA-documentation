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

colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")

## gather data ##
groupsPeaks <-"/path/to/peaks/"

cellFractions <- c("IRF8", "SOX")
df.list <- list()
for(cf in cellFractions){
  ## Load file with peaks grouped by overlap with samples
  load(paste0(groupsPeaks, cf,"_5e-2.narrowPeak.rdata"))
  df <-  as.data.frame(peak_set_group)
  df$pvals <- 10^(-df$pValue)
  df$qvals<- 10^(-df$qValue)
  df$fraction <- cf
  df.list[[cf]] <- df
}

peaks.df <- do.call(rbind,df.list)
peaks.df$peakID <- rownames(peaks.df)
rownames(peaks.df) <- NULL 
peaks.df$fraction <- ifelse(peaks.df$fraction == "IRF8", "IRF8+", "Sox10+")

# ########################## Fig 7.a ###########################################
# Violin plots of distribution of number of samples overlap with (y-axis)
# q-value threshold (x-axis). 
# coloured by cell-fraction (facets)

thresholds <- c(1e-100, 10^seq(-2, -100, length.out = 50), 0.05, 1)
thresholds <- sort(unique(thresholds))
order_levels.og <- as.character(head(thresholds, -1))

peaks.df$qValuegroup <- cut(
  peaks.df$qvals,
  breaks = thresholds,
  labels = order_levels.og,
  include.lowest = TRUE
)
peaks.df$qValuegroup[peaks.df$qvals < 1e-100] <- "1e-100"
peaks.df$qValuegroup <- factor(
  peaks.df$qValuegroup,
  levels = rev(levels(peaks.df$qValuegroup))
)
fig7a <- ggplot(peaks.df, aes(x = qValuegroup, y = n_samples, fill = fraction))+ 
  geom_violin(scale = "width")+
  scale_fill_manual(values = col_vals)+
  facet_wrap(~fraction, scales = "free_y")+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+
  labs(x = "Q value", y = "Number of samples overlap", title = "A") +
  pub_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45,hjust = 1))

# ########################## Fig 7.b ###########################################
# Line graph of number of peaks (y-axis)
# q-value threshold (x-axis) 
# with one line per “peak size group”.  
# coloured by cell-fraction (facets)

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

plot_df <- peaks.df %>% dplyr::count(qValuegroup, group, fraction)

fig7b <-ggplot(plot_df, aes(x = qValuegroup, y = n, group=group, color = group)) +
  geom_point() +
  geom_line(linewidth = 0.8)+
  facet_wrap(~fraction, scales = "free_y")+
  scale_y_continuous(
                     labels = scales::label_number(suffix = "K",scale = 1e-3))+
  labs(
    x = "Q value",
    y = "Number of peaks", title = "B", color = "Peak width\ngroup"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45,hjust = 1))

# ########################## Fig 7.c ###########################################
# Line graph of number of peaks (y-axis) 
# q-value threshold (x-axis) with
# one line per overlap with 0,1,2,3 samples etc.  
# plot 0 up to 10 individually, then 10-20, 20-30, or more than 30

peaks.df$n_samples_group <- cut(
  as.numeric(peaks.df$n_samples),
  breaks = c(0,1,2,3,4,5,6,7,8,9, 10, 20, 30, Inf),
  labels = c("0","1","2","3","4","5","6","7","8","9","10–20", "20–30", ">30"),
  right = TRUE,
  include.lowest = TRUE
)

table(peaks.df$qValuegroup)   # q-value bins (factor)
table(peaks.df$n_samples_group)     # number of samples per peak
peaks.df$qValuegroup <- factor(peaks.df$qValuegroup)
peaks.df$n_samples_group <- factor(peaks.df$n_samples_group)
tab <- table(peaks.df$qValuegroup, peaks.df$n_samples_group, peaks.df$fraction)
df_plot <- as.data.frame(tab)
colnames(df_plot) <- c("qValuegroup", "n_samples", "fraction", "count")

fig7c <-ggplot(df_plot, aes(x = qValuegroup, y = count, group=n_samples, color = n_samples)) +
  geom_point() +
  geom_line(linewidth = 1)+
  facet_wrap(~fraction, scales = "free_y")+
  labs(
    x = "Q value",
    y = "Number of peaks", title = "C", color = "Number of samples\noverlap"
  ) +scale_y_continuous(
    labels = scales::label_number(suffix = "K",scale = 1e-3))+
  pub_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## combined plot for Supplementary Figure 7 ##
design <- "AAAAA
           AAAAA
           BBBBB
           BBBBB
           CCCCC
           CCCCC"

combined_figSup7 <- fig7a + fig7b + fig7c +  plot_layout(design = design)

ggsave(
  filename = paste0("SupplementaryFig_7_dataset1.png"),
  plot = combined_figSup7,
  width = 24,
  height = 16
)

ggsave(
  filename = paste0("SupplementaryFig_7_dataset1.pdf"),
  plot = combined_figSup7,
  width = 24,
  height = 16
)
