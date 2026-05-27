#####################################################
##              Supplementary Figure 6             ##
##     Lack of reproducibility across samples      ##
#####################################################

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(reshape2)
library(patchwork)
library(tibble)

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

## Load file with peaks grouped by overlap with samples
groupsPeaks <-"/path/to/peaks/"
load(paste0(groupsPeaks,"NEUN_5e-2.subReads.rdata"))
names(allRes) <- c("tableCounts","samplesOverlap")
sampsOver <- allRes$samplesOverlap
tableCounts <- allRes$tableCounts
numPeaks <- sapply(tableCounts,length)

tableCounts.merged <- unlist(GRangesList(tableCounts), use.names = FALSE)
mcols(tableCounts.merged)$group <- rep(names(tableCounts), elementNROWS(GRangesList(tableCounts)))
names(tableCounts.merged) <- NULL
tableCounts.merged.df <- as.data.frame(tableCounts.merged)
tableCounts.merged.df$reads <- as.integer(str_extract(tableCounts.merged.df$group, "(?<=_)[0-9]+(?=readsAll)"))
colnames(tableCounts.merged.df) <- c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal",
                                     "pvalue", "qvalue", "score2", "n_samples", "peakFile","setReads")
tableCounts.merged.df$pvals <- 10^(-tableCounts.merged.df$pvalue)
tableCounts.merged.df$qvals<- 10^(-tableCounts.merged.df$qvalue)

# ########################## Fig Sup 6.a ###########################################
# Violin plot of number of samples each peak overlaps with (y-axis) 
# number of reads (x-axis) for default thresholds

fig6a <- ggplot(tableCounts.merged.df, aes(x = factor(setReads), y = n_samples))+ 
  geom_violin(scale = "width", fill = "#56B4E9")+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+scale_x_discrete(
    labels = function(x) paste0(x, "M")
  ) +
  labs(x = "Number of reads pairs", y = "Number of samples overlap", title = "A") +
  pub_theme +
  theme(
    legend.position = "none")

# ########################## Fig Sup 6.b ###########################################
# Line graph of mean number of samples each peak overlaps with (y-axis) 
# number of reads (x-axis) for default thresholds
# with separate lines for different groups of peak sizes (i.e. very narrow, narrow,…)

tableCounts.merged.df$widthgroup <- cut(
  tableCounts.merged.df$width,
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
tableCounts.merged.df$setReadsfact <- factor(tableCounts.merged.df$setReads)
plot_df <- tableCounts.merged.df %>%
  group_by(setReadsfact, widthgroup) %>%
  summarise(mean_n_sample = mean(n_samples, na.rm = TRUE), .groups = "drop")

fig6b <-ggplot(plot_df, aes(x = factor(setReadsfact), y = mean_n_sample, group=widthgroup,
                         color = widthgroup )) +
  geom_point() +
  geom_line(linewidth = 0.8)+
  scale_x_discrete(
    labels = function(x) paste0(x, "M")
  )+
  labs(x = "Number of reads pairs", y = "Mean samples overlap",
       title = "B", color = "Peak width group"
  ) +pub_theme +
  theme(
    legend.position = "right"
  )

# ########################## Fig Sup 6.c ###########################################
# Scatterplot of number of peaks (x-axis) 
# proportion of peaks overlapping merged bam peakset (one point per sample) 


sampsNames <- colnames(sampsOver[[1]])

df <- imap_dfr(sampsOver, function(df, peak_set_name) {
  
  df <- as.data.frame(df)
  
  n_row <- grep("n_total", rownames(df))
  p_row <- grep("pct_overlap", rownames(df))

  n_total <-as.integer(t(df[n_row, ]))
  pct_overlap <-as.integer(t(df[p_row, ]))
  
  
 df <-  data.frame(peak_set = peak_set_name,sample = colnames(df),n_total = n_total,pct_overlap = pct_overlap)
 rownames(df) <- NULL
  df
})

df$setReads <- stringr::str_extract(df$peak_set, "(?<=_)[0-9]+(?=readsAll)")

fig6c <-ggplot(df, aes(x = round(n_total), y = pct_overlap, colour = setReads)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0,300000,by = 50000),
                     labels = scales::label_number(suffix = "K",scale = 1e-3))+
  labs(x = "Number of peaks", y = "% peaks overlap peak set",
       title = "C", color = "Number of reads\nmerged bam (M)"
  ) +pub_theme +
  theme(legend.text = element_text(size = 11),
        legend.position = "bottom"
  )+guides(colour = guide_legend(override.aes = list(size=6)))

# ########################## Fig Sup 6.d ###########################################
# Scatterplot of number of reads (x-axis)
# proportion of peaks overlapping merged bam peakset (one point per sample) 

## read file with number of reads for each sample
samplesReads <- read.table(paste0(groupsPeaks,"fragment_counts_all_NEUN_reads_BDR.txt"), sep = "\t")
samplesReads$sample <- sub("dir_string", "", samplesReads$V1)
samplesReads$sample <- gsub(".filt.nodup.bam", "", samplesReads$sample)

df$numReads <- samplesReads[match(df$sample,samplesReads$sample),]$V2
df <- df[complete.cases(df$numReads),]

fig6d <-ggplot(df, aes(x = numReads, y = pct_overlap, color = setReads)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0,400000000,by = 100000000),
                     labels = scales::label_number(suffix = "M",scale = 1e-6))+
  labs(x = "Number of reads pairs", y = "% peaks overlap peak set",
       title = "D", color = "Number of reads\nmerged bam (M)"
  ) +pub_theme +
  theme(legend.text = element_text(size = 11),
        legend.position = "bottom"
  )+guides(colour = guide_legend(override.aes = list(size=6)))

legend_shared_bottom <- cowplot::get_legend(
  fig6d +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",legend.justification = "center"
    )
)

legend_shared_upper <- cowplot::get_legend(
  fig6b +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
)

fig6d_clean <- fig6d + theme(legend.position = "none")
fig6c_clean <- fig6c + theme(legend.position = "none")
fig6b_clean <- fig6b + theme(legend.position = "none")


## Collate plots for combined Supplementary Figure 6 ##
design <- "AAABBB
           AAABBB
           CCCDDD
           CCCDDD
           ...L..
           ...L.."

# Put legend inside a centered drawing canvas
legend_centered <- ggdraw() +
  draw_grob(legend_shared_bottom, x = 0.5, hjust = 0.5)

combined_fig6 <-
  fig6a + fig6b + fig6c_clean +fig6d_clean + legend_centered +
  plot_layout(
    heights = c(1,1,0.15),
    nrow = 3
  )
ggsave(
  filename = "SupplementaryFig_6_dataset1_NEUN.pdf",
  plot = combined_fig6,
  width = 16,
  height = 12
)

ggsave(
  filename = "SupplementaryFig_6_dataset1_NEUN.png",
  plot = combined_fig6,
  width = 16,
  height = 12
)
