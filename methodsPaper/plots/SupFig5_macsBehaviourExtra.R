#######################################
##       Figure Supplementary 5     ##
#######################################

library(ggplot2)
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

col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")
peakFiles<-list.files(groupsPeaks, pattern="SOX_1e-2.*\\.narrowPeak\\.filt|SOX_5e-2.*\\.narrowPeak\\.filt|SOX_1e-5.*\\.narrowPeak\\.filt")
print("Peaks files: ")
print(peakFiles)

peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(file.path(groupsPeaks, each)), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
}
numPeaks <- sapply(peaks, length)
widthPeaks <- sapply(peaks, width)
meanWidthPeaks<- sapply(widthPeaks, mean)

df <- tibble(file = names(numPeaks)) %>%
  mutate(
    num = numPeaks,
    cutoff = str_extract(file, "[0-9]+e-[0-9]+"),
    reads    = as.integer(str_extract(file, "(?<=_)[0-9]+(?=readsAll)"))
  )

df$cutoff <-fct_recode(df$cutoff, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")

# ########################## Fig 4.a ###########################################
# Line graph of number of peaks (y-axis)
# number of reads (x-axis) 
# coloured by default thresholds

p1 <- ggplot(df, aes(x = reads, y = num, color = cutoff)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_color_manual(values = col_vals) +
  scale_x_continuous(
    breaks = sort(unique(df$reads)),
    labels = scales::label_number( suffix = "M")
  ) +
  scale_y_continuous(
    labels = scales::label_number( suffix = "K",scale = 1e-3)
  )+
  labs(
    x = "Number of reads pairs",
    y = "Number of peaks",
    title = "A",
    color = "Q value threshold"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig 4.b ###########################################
# Line graph of number of peaks (y-axis)
# number of reads (x-axis) for default thresholds
# with separate lines for different groups of peak sizes

peaksAnno.all.df.space$group <- cut(
  peaksAnno.all.df.space$width,
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

plot_df <- peaksAnno.all.df.space %>% count(set_num, group)
p2 <-ggplot(plot_df, aes(x = set_num, y = n, group=group, color = group)) +
  geom_point() +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(
    breaks = sort(unique(plot_df$set_num)),
    labels = scales::label_number( suffix = "M")
  ) +
  scale_y_continuous(
    labels = scales::label_number( suffix = "M",scale = 1e-4)
  )+
  labs(
    x = "Number of reads",
    y = "Number of peaks", title = "B", color = "Peak width\ngroup"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom")

# ########################## Fig 4.c ###########################################
# Box plot/violin plots of width of peaks (y-axis)
# number of reads (x-axis) 
# for peak set called with 5e-2

# subset to only peak sets with 5e-2 q-value threshold
peak_sets_sub <- peaks[grepl("5e-2", names(peaks))]
widthPeaks_vec <- sapply(peak_sets_sub, width)
numPeaks.sub <- sapply(peak_sets_sub, length)
widthPeaks_df <- data.frame(
  set = rep(names(numPeaks.sub), times = numPeaks.sub),
  width = unlist(widthPeaks_vec)
)

widthPeaks_df$set_num <- as.integer(str_extract(widthPeaks_df$set, "(?<=_)[0-9]+(?=readsAll)"))

p3 <-ggplot(widthPeaks_df, aes(x = factor(set_num), y = log10(width)))+ 
  geom_violin(scale = "width", fill = "#56B4E9") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white") +
  labs(x = "Number of reads pairs", y = "base pairs (scale log10)", title = "B") +
  pub_theme +
  theme(
    legend.position = "none"
  ) + scale_x_discrete(
    labels = function(x) paste0(x, "M")
  )

# ########################## Fig 4.d ###########################################
# Box plot/violin plots of gaps between peaks (y-axis) 
# number of reads (x-axis) 
# for peak set called with 5e-2

spacing_prev <- lapply(peaks, function(x) {
  hits <- distanceToNearest(x, ignore.strand = TRUE)
  dist <- rep(NA_integer_, length(x))
  dist[queryHits(hits)] <- mcols(hits)$distance
  
  dist
})

peaksAnno.all.df.space <- bind_rows(peaksAnno.all.df)

peaksAnno.all.df.space$set_num <- as.integer(str_extract(peaksAnno.all.df.space$set, "(?<=_)[0-9]+(?=readsAll)"))

p4 <- ggplot(peaksAnno.all.df.space, aes(x = factor(set_num), y = log10(spacing)))+ 
  geom_violin(scale = "width", fill = "#56B4E9")+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")  +
  labs(x = "Number of reads pairs", y = "base pairs to nearest peak (scale log10)",title = "C") +
  pub_theme +
  theme(
    legend.position = "none"
  ) + scale_x_discrete(
    labels = function(x) paste0(x, "M")
  )

## combined figure with all plots
design <- "AAABBB
           AAABBB
           DDDCCC
           DDDCCC"
combined_figSup4 <- p1 + p2 + p3+ p4 +  plot_layout(design = design)

ggsave(
  filename = "SupplementaryFig_4_dataset1.pdf",
  plot = combined_figSup4,
  width = 15,
  height = 12
)
