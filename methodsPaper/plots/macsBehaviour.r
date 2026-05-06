library(ggplot2)
library(GenomicRanges)

peakFiles<-list.files(groupsPeaks, pattern="NEUN_5e-2.*\\.narrowPeak\\.filt")
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
    cutoff = str_extract(file, "[0-9]+e-[0-9]+"),
    reads    = as.integer(str_extract(file, "(?<=_)[0-9]+(?=readsAll)"))
  )

p1<-ggplot(df, aes(x = reads, y = num / 1000)) +
  scale_x_continuous(breaks = sort(unique(df$reads))) +  # guarantees numeric ordering) 
  geom_point(size = 3, color = '#FF6A6A')+
  geom_line(size = 1, color = '#FF6A6A')+
  labs(x = "Number of reads pairs (x10^8)",
       y = "Number of peaks (x1000)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 11),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold")
  )

widthPeaks_vec <- sapply(peaks, width)
widthPeaks_df <- data.frame(
  set = rep(names(numPeaks), times = numPeaks),
  width = unlist(widthPeaks_vec)
)

widthPeaks_df$set_num <- as.integer(str_extract(widthPeaks_df$set, "(?<=_)[0-9]+(?=readsAll)"))

p2 <-ggplot(widthPeaks_df, aes(x = factor(set_num), y = log10(width)))+ 
  geom_violin(scale = "width", fill = '#FF6A6A') +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white") +
  labs(x = "Number of reads pairs (x10^8)", y = "Width (log10bp)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 11),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold")
  )

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

p3 <- ggplot(peaksAnno.all.df.space, aes(x = factor(set_num), y = log10(spacing)))+ 
  geom_violin(scale = "width", fill = '#FF6A6A')+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")  +
  labs(x = "Number of reads pairs (x10^8)", y = "Distance to nearest peak (log10bp)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 11),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold")
  )

design <- "#AAAA#
           #AAAA#
           BBBCCC
           BBBCCC"
combined_plot <- free(p1) + p2 + p3 +  plot_layout(axes = "collect",design = design)

# Show combined plot
combined_plot
