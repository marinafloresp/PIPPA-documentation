library(tidyr)
library(dplyr)
library(GenomicRanges)

## Universal theme for plotting
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )

## specify cell-fraction
cf <- "NEUN"

## Load file with peaks grouped by overlap with samples
load(paste0(cf,"_5e-2.narrowPeak.rdata"))
head(peak_set_group)

peaks.df <- as.data.frame(peak_set_group)
peaks.df$pvals <- 10^(-peaks.df$pValue)
peaks.df$qvals<- 10^(-peaks.df$qValue)

# ########################## Fig 3b.a ###########################################
# Violin plots of distribution of number of samples overlap with (y-axis)
# q-value threshold (x-axis). 

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

p1 <- ggplot(peaks.df, aes(x = qValuegroup, y = n_samples))+ 
  geom_violin(scale = "width", fill = '#FF6A6A')+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+
  labs(x = "Q-value", y = "Number of samples overlap", title = "A") +
  pub_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45,hjust = 1))

# ########################## Fig 3b.b ###########################################
# Line graph of number of peaks (y-axis)
# q-value threshold (x-axis) 
# with one line per “peak size group”.  

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

plot_df <- peaks.df %>% count(qValuegroup, group)

#line-graph version of plot
p2 <-ggplot(plot_df, aes(x = qValuegroup, y = n, group=group, color = group)) +
  geom_point() +
  geom_line(linewidth = 0.8)+
  labs(
    x = "Q value",
    y = "Number of peaks", title = "B", color = "Peak width\ngroup"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom",
     axis.text.x = element_text(angle = 45,hjust = 1))


# ########################## Fig 3b.c ###########################################
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
tab <- table(peaks.df$qValuegroup, peaks.df$n_samples_group)
df_plot <- as.data.frame(tab)
colnames(df_plot) <- c("qValuegroup", "n_samples", "count")

#line-graph version of plot
p3 <-ggplot(df_plot, aes(x = qValuegroup, y = count, group=n_samples, color = n_samples)) +
  geom_point() +
  geom_line(size = 1)+
  labs(
    x = "Q value",
    y = "Number of peaks", title = "C", color = "Number of samples\noverlap"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## Combine plots into single figure ##
combined_fig3b <-
  (p1/p2 / p3) + 
  plot_layout(heights = c(1, 1, 1))
combined_fig3b

## Save figure ##
ggsave(
  filename = "Figure_3.b_dataset1_NEUN.pdf",
  plot = combined_fig3b,
  width = 12,
  height = 16
)
ggsave(
  filename = "Figure_3.b_dataset1_NEUN.png",
  plot = combined_fig3b,
  width = 12,
  height = 16
)

# stacked barplot version #
p2 <- ggplot(plot_df, aes(x = qValuegroup, y = n, group=group, fill = group)) +
  geom_bar(stat = "identity") + 
  labs(
    x = "Q value",
    y = "Number of peaks", title = "B", fill = "Peak width group"
  ) +pub_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )
p3 <- ggplot(df_plot, aes(x = qValuegroup, y = count, fill = n_samples)) +
  geom_bar(stat = "identity") + 
  labs(
    x = "Q value",
    y = "Number of peaks", title = "C", fill = "Number of samples\noverlap"
  ) +pub_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## Combine plots into single figure ##
combined_fig3b <-
  (p1/p2 / p3) + 
  plot_layout(heights = c(1, 1, 1))
combined_fig3b

## Save figure ##
ggsave(
  filename = "Figure_3.b_bars_dataset1_NEUN.pdf",
  plot = combined_fig3b,
  width = 12,
  height = 16
)
ggsave(
  filename = "Figure_3.b_bars_dataset1_NEUN.png",
  plot = combined_fig3b,
  width = 12,
  height = 16
)
