########################################
##       Supplementary Fig 11,12      ##
########################################

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
    plot.title = element_text(face = "bold", size = 24),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")
colorBlindAlt <- c(
  "#332288",  # dark blue/purple
  "#88CCEE",  # light blue
  "#44AA99",  # teal
  "#117733",  # green
  "#DDCC77",  # sand
  "#CC6677",  # muted red
  "#AA4499",  # purple
  "#882255"   # dark magenta
)

## Load file with peaks grouped by overlap with samples
load("./data/downReads/NEUN_5e-2.subReads.rdata")
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

# ########################## Fig Sup 11.b ###########################################
# Stacked bar plot graph of % from total peaks (y-axis) 
# number of reads (x-axis) for default thresholds
# with separate groups based on number of samples peaks overlap with

plot_df <- tableCounts.merged.df %>%
  dplyr::mutate(
    overlap_group = case_when(
      n_samples == 0 ~ "0",
      n_samples == 1 ~ "1",
      n_samples == 2 ~ "2",
      n_samples >= 3  & n_samples <= 6  ~ "3-6",
      n_samples >= 7  & n_samples <= 10 ~ "7-10",
      n_samples >= 11 & n_samples <= 20 ~ "11-20",
      n_samples >= 21 & n_samples <= 30 ~ "21-30",
      n_samples > 30 ~ ">30"
    ),
    overlap_group = factor(
      overlap_group,
      levels = c("0", "1", "2", "3-6", "7-10",
                 "11-20", "21-30", ">30")
    )
  ) %>%
  dplyr::count(setReads, overlap_group)%>%
  dplyr::group_by(setReads) %>%
  dplyr::mutate(percentage = n / sum(n) * 100) %>%
  dplyr::ungroup()

plot_df[plot_df$setReads == 100,]

plot_df <- plot_df %>%
  dplyr::mutate(
    overlap_group = factor(
      overlap_group,
      levels = c(
        "0", "1", "2", "3-6",
        "7-10", "11-20", "21-30", ">30"
      )
    )
  )
p11b <- ggplot(plot_df, aes(
    x = factor(setReads),
    y = percentage,
    fill = overlap_group
  )) +
  scale_fill_manual(values = colorBlindGrey8)+
  geom_col(width = 0.95,position = position_stack(reverse = TRUE)) +
  pub_theme + 
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "% of total peaks",
    fill =  "Number of samples\noverlap", title = "B"
  ) + theme(legend.position = "none")


# ########################## Fig Sup 11.c ###########################################
# Stacked bar plot graph of number of peaks (y-axis) 
# number of reads (x-axis) for default thresholds
# with separate groups based on number of samples peaks overlap with

plot_df <- tableCounts.merged.df %>%
  dplyr::mutate(
    overlap_group = case_when(
      n_samples == 0 ~ "0",
      n_samples == 1 ~ "1",
      n_samples == 2 ~ "2",
      n_samples >= 3  & n_samples <= 6  ~ "3-6",
      n_samples >= 7  & n_samples <= 10 ~ "7-10",
      n_samples >= 11 & n_samples <= 20 ~ "11-20",
      n_samples >= 21 & n_samples <= 30 ~ "21-30",
      n_samples > 30 ~ ">30"
    ),
    overlap_group = factor(
      overlap_group,
      levels = c("0", "1", "2", "3-6", "7-10",
                 "11-20", "21-30", ">30")
    )
  ) %>%
  dplyr::count(setReads, overlap_group)%>%
  dplyr::group_by(setReads) %>%
  dplyr::mutate(partNum = n ) %>%
dplyr::ungroup()

plot_df[plot_df$setReads == 100,]

plot_df <- plot_df %>%
  dplyr::mutate(
    overlap_group = factor(
      overlap_group,
      levels = c(
        "0", "1", "2", "3-6",
        "7-10", "11-20", "21-30", ">30"
      )
    )
  )


p11c <- ggplot(plot_df, aes(
  x = factor(setReads),
  y = partNum,
  fill = overlap_group
)) +scale_fill_manual(values = colorBlindGrey8)+
  geom_col(width = 0.95,position = position_stack(reverse = TRUE)) +
  pub_theme + 
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Number of peaks",
    fill =  "Number of samples\noverlap", title = "C"
  )+ theme(legend.position = "bottom")



# ########################## Fig Sup 11.a ###########################################
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
  dplyr::group_by(setReadsfact, widthgroup) %>%
  dplyr::summarise(mean_n_sample = mean(n_samples, na.rm = TRUE), .groups = "drop")

##stats for paper
plot_df[plot_df$setReadsfact == 100,]

#line-graph version of plot
p11a <-ggplot(plot_df, aes(x = factor(setReadsfact), y = mean_n_sample, group=widthgroup,
                            color = widthgroup )) +
  geom_point() +
  geom_line(linewidth = 1)+
  scale_x_discrete(
    labels = function(x) paste0(x)
  )+scale_color_manual(values = colorBlindAlt)+
  labs(x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), y = "Mean samples overlap",
       title = "A", color = "Peak width group"
  ) +pub_theme +
  theme(
    legend.position = "bottom"
  )
p11a

legend_A <- get_legend(
  p11a +
    theme(legend.position = "bottom")
)

# Extract B/C legend
legend_BC <- get_legend(
  p11b +
    theme(legend.position = "bottom")
)

# Remove legends from plots
p11a_noleg <- p11a +
  theme(legend.position = "none")

p11b_noleg <- p11b +
  theme(legend.position = "none")

p11c_noleg <- p11c +
  theme(legend.position = "none")

# Main figure
main_plots <- p11a_noleg | p11b_noleg | p11c_noleg

# Legends side-by-side
legends <- wrap_elements(legend_A) |
  wrap_elements(legend_BC)

# Combine
supFig11 <- main_plots /
  legends +
  plot_layout(
    heights = c(1, 0.15)
  )
ggsave(
  filename = "SupplementaryFig_11_D1_NEUN.pdf",
  plot = supFig11,
  width = 22,
  height = 6
)

ggsave(
  filename = "SupplementaryFig_11_D1_NEUN.png",
  plot = supFig11,
  width = 16,
  height = 6)

# ########################## Fig Sup 12.a ###########################################
# Violin plot of number of samples each peak overlaps with (y-axis) 
# number of reads (x-axis) for default thresholds

p12a <- ggplot(tableCounts.merged.df, aes(x = factor(setReads), y = n_samples))+ 
  geom_violin(scale = "width", fill = "#56B4E9")+
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")+
  scale_x_discrete(
    labels = function(x) paste0(x)
  ) +
  labs(x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), y = "Number of samples overlap", title = "A") +
  pub_theme +
  theme(
    legend.position = "none")


load("./methodsPaper/data/downReads/NEUN_5e-2.subReads1overlap.rdata")
allRes.unique <- allRes$uniqueSummary

plot_df <- data.frame(
  sample = allRes.unique$sample,
  sample_peaks = as.numeric(allRes.unique$sample_peaks),
  unique_group_peaks = as.numeric(allRes.unique$unique_group_peaks),
  group = allRes.unique$group
)

samplesReads <- read.table("./methodsPaper/data/reads_counts_passQC_NEUN_reads.txt", sep = "\t")
samplesReads$sample <- gsub(".filt.nodup.bam", "", samplesReads$sample)
samplesReads <- samplesReads[,c(2,3)]
plot_df$numReads <- samplesReads[match(plot_df$sample,samplesReads$sample),]$V2

# Fraction of each sample's peaks that are uniquely represented
plot_df$unique_fraction <- 
  plot_df$unique_group_peaks / plot_df$sample_peaks

plot_df.sub <- plot_df %>% filter(group == "50")

# ########################## Fig Sup 12.b ###########################################
# Scatter plot of number of peaks that only overlap one sample (y-axis)
# by the number of reads for that sample they overlap with only (x-axis) 

p12b <- ggplot(plot_df.sub, aes(
  x = numReads,
  y = unique_group_peaks)) +
  geom_point(
    size = 3,color = "#56B4E9")+
  scale_x_continuous(breaks = seq(0,max(plot_df.sub$numReads),by = 50000000),
                     labels = scales::label_number(suffix = "",scale = 1e-6))+
  pub_theme +
  labs(
    x = expression("Number of reads per sample (" ~ x ~ 10^6 ~ ")"),
    y = "Peaks overlap one sample", title = "B"
  ) 

# ########################## Fig Sup 12.c ###########################################
# Scatter plot of number of peaks per sample (y-axis)
# by number of peaks that only overlap one sample (x-axis) 

p12c <- ggplot(plot_df.sub, aes(
  x = sample_peaks ,
  y = unique_group_peaks)) +
  geom_point( size = 3,color = "#56B4E9")+
  scale_x_continuous(breaks = seq(0,max(plot_df.sub$sample_peaks),by = 50000),
                     labels = scales::label_number(suffix = "",scale = 1e-3))+
  pub_theme + theme(legend.position = "right")+
  labs(
    x = expression("Peaks per sample (" ~ x ~ 10^3 ~ ")"),
    y = "Peaks overlap one sample", title = "C", colour = ""
  )

supFig12 <- (p12a + p12b +p12c)+
  plot_layout(
    heights = c(1,1),
    nrow = 1
  )
ggsave(
  filename = "SupplementaryFig_12_D1_NEUN.pdf",
  plot = supFig12,
  width = 22,
  height = 6
)

ggsave(
  filename = "SupplementaryFig_12_D1_NEUN.png",
  plot = supFig12,
  width = 16,
  height = 6)
