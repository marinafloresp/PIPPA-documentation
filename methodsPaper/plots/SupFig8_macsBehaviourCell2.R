#####################################
##       Supplementary Fig 8       ##
#####################################

library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(patchwork)
library(cowplot)
options(scipen = 5)

## settings for plotting ##

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
cf <- "SOX"
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")
col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")

######################
## Gather peak sets ##
######################

groupsPeaks <-"./methodsPaper/data/downReads/"
peakFiles<-list.files(groupsPeaks, pattern=paste0(cf,"_1e-2.*\\.narrowPeak\\.filt|",cf,"_5e-2.*\\.narrowPeak\\.filt|",cf,"_1e-5.*\\.narrowPeak\\.filt"))
files_readsAll <- peakFiles[grepl("_\\d+readsAll\\.narrowPeak\\.filt", peakFiles)]
files_readsAll <- files_readsAll[grepl(cf, files_readsAll, ignore.case = TRUE)]

print("Peaks files: ")
print(files_readsAll)

peaks<-list()
for(each in files_readsAll){
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

# ########################## Fig Sup 8.a ###########################################
# Line plot of number of peaks (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold

p8a <- ggplot(df, aes(x = reads, y = num, color = cutoff)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_vals) +
  scale_x_continuous(
    breaks = sort(unique(df$reads))
  ) +
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Number of peaks",
    title = "A",
    color = "Q value threshold"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom", plot.title = element_text(size = 24)
  )


# ########################## Fig Sup 8.c ###########################################
# Violin plot of width of peaks (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold ( default)

peak_sets_sub <- peaks[grepl("5e-2", names(peaks))]
widthPeaks_vec <- sapply(peak_sets_sub, width)
numPeaks.sub <- sapply(peak_sets_sub, length)
widthPeaks_df <- data.frame(
  set = rep(names(numPeaks.sub), times = numPeaks.sub),
  width = unlist(widthPeaks_vec)
)

widthPeaks_df$set_num <- as.integer(str_extract(widthPeaks_df$set, "(?<=_)[0-9]+(?=readsAll)"))
sum(widthPeaks_df[widthPeaks_df$set_num == 100,]$width)

p8c <-ggplot(widthPeaks_df, aes(x = factor(set_num), y = log10(width)))+ 
  geom_violin(scale = "width", fill = "#56B4E9") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white") +
  labs(x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), y = "Peak width (bp, log10 scale)", title = "C") +
  pub_theme +
  theme(
    legend.position = "none", plot.title = element_text(size = 24)
  ) 

# ########################## Fig Sup 8.d ###########################################
# Violin plot of space between peaks (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold ( default)

spacing_prev <- lapply(peak_sets_sub, function(x) {
  hits <- distanceToNearest(x, ignore.strand = TRUE)
  dist <- rep(NA_integer_, length(x))
  dist[queryHits(hits)] <- mcols(hits)$distance
  
  dist
})


peaksAnno.all.df <- lapply(seq_along(peak_sets_sub), function(i) {
  df <- as.data.frame(peak_sets_sub[[i]])
  df$spacing <- spacing_prev[[i]]
  df$set <- names(peak_sets_sub)[i]
  df
})

peaksAnno.all.df.space <- bind_rows(peaksAnno.all.df)

peaksAnno.all.df.space$set_num <- as.integer(str_extract(peaksAnno.all.df.space$set, 
                                                         "(?<=_)[0-9]+(?=readsAll)"))

p8d <- ggplot(peaksAnno.all.df.space, aes(x = factor(set_num), y = log10(spacing)))+ 
  geom_violin(scale = "width", fill = "#56B4E9")+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4,fill = 'black') +
  stat_summary(fun = median, geom = "point", color = "white")  +
  labs(x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"), 
  y = "Gap between peaks (bp, log10 scale)",title = "D") +
  pub_theme +
  theme(
    legend.position = "none", plot.title = element_text(size = 24)
  ) 
p8d

# ########################## Fig Sup 8.b ###########################################
# Line plot of number of peaks (y-axis)
# by read pairs (x-axis) 
# coloured by peak width group

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

plot_df <- peaksAnno.all.df.space %>% dplyr::count(set_num, group)

p8b <-ggplot(plot_df, aes(x = set_num, y = n, group=group, color = group)) +
  geom_point() +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = sort(unique(plot_df$set_num))
    
  ) + scale_color_manual(values = colorBlindGrey8)+
  labs(
    x =expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Number of peaks", title = "B", color = "Peak width\ngroup"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom", plot.title = element_text(size = 24))

design <- "AAADDD
           AAADDD
           BBBCCC
           BBBCCC"
figSup8 <- p8a + p8b + p8c+ p8d +  plot_layout(design = design)

ggsave(
  filename = "SupplementaryFig_8_D1_SOX.pdf",
  plot = figSup8,
  width = 15,
  height = 12
)

ggsave(
  filename = "SupplementaryFig_8_D1_SOX.png",
  plot = figSup8,
  width = 15,
  height = 12
)
