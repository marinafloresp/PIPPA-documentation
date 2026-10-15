#########################
##       Figure 2      ##
#########################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)

groupsPeaks <- "/dir/to/peak/sets"

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

cf <- "NEUN"
peakFiles<-list.files(groupsPeaks, pattern=paste0(cf,"_1e-2.*\\.narrowPeak\\.filt|",cf,"_5e-2.*\\.narrowPeak\\.filt|",cf,"_1e-5.*\\.narrowPeak\\.filt"))
files_readsAll <- peakFiles[grepl("_\\d+readsAll\\.narrowPeak\\.filt", peakFiles)]
files_readsAll <- files_readsAll[grepl(cf, files_readsAll, ignore.case = TRUE)]

print("Peaks files: ")
print(files_readsAll)

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

# ########################## Fig 2.a ###########################################
# Line graph of number of peaks (y-axis)
# number of reads (x-axis) 
# coloured by default thresholds 

p2a <- ggplot(df, aes(x = reads, y = num, color = cutoff)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_vals) +
  scale_x_continuous(
    breaks = sort(unique(df$reads))
  )   +
labs(
  x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    #x = "Read pairs (Million)",
    y = "Number of peaks",
    title = "A",
    color = "Q value threshold"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom", plot.title = element_text(size = 24)
  )

# ########################## Fig 2.b ###########################################
# Box plot/violin plots of width of peaks (y-axis) 
# number of reads (x-axis) 
# only 5e-2 peaks for all read subsets

# subset to only peak sets with 5e-2 q-value threshold
peak_sets_sub <- peaks[grepl("5e-2", names(peaks))]
widthPeaks_vec <- sapply(peak_sets_sub, width)
numPeaks.sub <- sapply(peak_sets_sub, length)
widthPeaks_df <- data.frame(
  set = rep(names(numPeaks.sub), times = numPeaks.sub),
  width = unlist(widthPeaks_vec)
)

widthPeaks_df$set_num <- as.integer(str_extract(widthPeaks_df$set, "(?<=_)[0-9]+(?=readsAll)"))

p2b <-ggplot(widthPeaks_df, aes(x = factor(set_num), y = width))+ 
  geom_violin(scale = "width", fill = "#56B4E9") +
  scale_y_log10()+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white") +
  labs( x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), y = expression("Peak width (bp; log"[10]*" scale)"),
  title = "B") +
  pub_theme +
  theme(
    legend.position = "none", plot.title = element_text(size = 24)
  ) 

# ########################## Fig 2.c ###########################################
# Box plot/violin plots of gaps between peaks (y-axis) 
# number of reads (x-axis) 
# only 5e-2 peaks for all read subsets

spacing_prev <- lapply(peaks, function(x) {
  hits <- distanceToNearest(x, ignore.strand = TRUE)
  dist <- rep(NA_integer_, length(x))
  dist[queryHits(hits)] <- mcols(hits)$distance
  
  dist
})

peaksAnno.all.df <- lapply(seq_along(peaks), function(i) {
  df <- as.data.frame(peaks[[i]])
  df$spacing <- spacing_prev[[i]]
  df$set <- names(peaks)[i]
  df
})

peaksAnno.all.df.space <- bind_rows(peaksAnno.all.df)

peaksAnno.all.df.space$set_num <- as.integer(str_extract(peaksAnno.all.df.space$set, "(?<=_)[0-9]+(?=readsAll)"))

p2c <- ggplot(peaksAnno.all.df.space, aes(x = factor(set_num), y = spacing))+ 
  geom_violin(scale = "width", fill = "#56B4E9")+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  scale_y_log10()+
  stat_summary(fun = median, geom = "point", color = "white")  +
  labs(x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), 
       y = expression("Gap between peaks (bp; log"[10]*" scale)"),title = "C") +
  pub_theme +
  theme(
    legend.position = "none", plot.title = element_text(size = 24)
  ) 

design <- "AAAAAA
           AAAAAA
           BBBCCC
           BBBCCC"

figure2 <- free(p2a) + p2b + p2c +  plot_layout(axes = "collect",design = design)

ggsave(
  filename = "Figure_2_D1_NEUN.pdf",
  plot = figure2,
  width = 12,
  height = 12
)

ggsave(
  filename = "Figure_2_D1_NEUN.png",
  plot = figure2,
  width = 12,
  height = 12
)
