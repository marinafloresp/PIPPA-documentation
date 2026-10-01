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

groupsPeaks <-"./methodsPaper/data/downReads/"
col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")

######################
## Gather peak sets ##
######################

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
