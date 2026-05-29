#####################################################
##       Supplementary Figure 5                    ##
## Genomic distribution of “inflated” peak sets at ## 
## default thresholds ( for all cell types).       ##
#####################################################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ggplot2)
library(tidyverse) 

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

groupsPeaks <-"/path/to/peaks/"
col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")

## Find all peak sets called in NeuN+, Sox10+ and IRF8+ with different q-value thresholds and number of random reads
peakFiles<-list.files(groupsPeaks, pattern="_1e-2.*\\.narrowPeak\\.filt|_5e-2.*\\.narrowPeak\\.filt|_1e-5.*\\.narrowPeak\\.filt")
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

peakAnno.all <- lapply(peaks, annotatePeak,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db")
names(peakAnno.all) <- names(peaks)

calculatePeaksGenes <- function(peakAnno, setName){
  
  peaks_per_gene <- peakAnno %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::count(SYMBOL, name = "peak_count")%>%
    dplyr::mutate(peak_set = setName)
  
  list(
    summary = tibble::tibble(
      name = setName,
      n_genes = nrow(peaks_per_gene),
      mean_peaks_per_gene = mean(peaks_per_gene$peak_count)
    ),
    peaks_per_gene = peaks_per_gene
  )
}

genesReadsSets <- Map(
  calculatePeaksGenes,
  peakAnno.all,
  names(peakAnno.all)
)
all_summaries <- lapply(genesReadsSets, `[[`, "summary")


genesReadsSets.all <- bind_rows(all_summaries)
genesReadsSets.all$num_reads <- as.integer(str_extract(genesReadsSets.all$name, "(?<=_)[0-9]+(?=readsAll)"))
genesReadsSets.all$thres <- str_extract(genesReadsSets.all$name, "[0-9]+e-[0-9]+")
genesReadsSets.all$fraction <-str_extract(genesReadsSets.all$name, "^[^_]+")
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "NEUN"] <- "NeuN+"
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "SOX"] <- "Sox10+"
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "IRF8"] <- "IRF8+"
genesReadsSets.all$thres <-fct_recode(genesReadsSets.all$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")


## ########################## Sup Fig 5.a ###########################################
#Line graph of number of genes annotated to any peak (y-axis) 
# number of reads (x-axis) for default thresholds#
# coloured by threshold
# linetype by cell-fraction

fig5a <- ggplot(as.data.frame(genesReadsSets.all), aes(x = num_reads, y = n_genes, color = thres, linetype = fraction)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSets.all$num_reads)),
    labels = scales::label_number( suffix = "M")
  ) +
  scale_linetype_manual(values = c("NeuN+" = "solid", "Sox10+" = "dotted", "IRF8+" = "dashed")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Number of reads pairs",
    y = "Number of annotated genes",
    title = "A",
    color = "Q value threshold", linetype = "Cell-fraction"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )


## ########################## Sup Fig 5.b ###########################################
# Violin plot of the number of peaks annotated to each gene (y-axis)
# number of reads (x-axis) for default thresholds. 
# column panels for thresholds
# row panels for cell-fractions

all_peaksGenes <- lapply(genesReadsSets, `[[`, "peaks_per_gene")
genesReadsSets.all <- bind_rows(all_peaksGenes)
genesReadsSets.all$num_reads <- as.integer(str_extract(genesReadsSets.all$peak_set, "(?<=_)[0-9]+(?=readsAll)"))
genesReadsSets.all$thres <- str_extract(genesReadsSets.all$peak_set, "[0-9]+e-[0-9]+")
genesReadsSets.all$thres <-fct_recode(genesReadsSets.all$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")

genesReadsSets.all$thres <- factor(
  genesReadsSets.all$thres,
  levels = c("0.05" ,"0.01" , "0.00001")
)
genesReadsSets.all$fraction <-str_extract(genesReadsSets.all$peak_set, "^[^_]+")
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "NEUN"] <- "NeuN+"
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "SOX"] <- "Sox10+"
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "IRF8"] <- "IRF8+"

fig5b2 <- ggplot(
  genesReadsSets.all,
  aes(x = factor(num_reads), y = peak_count, fill = thres)
) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    alpha = 0.3,
    fill = "black"
  ) +scale_y_log10()+
  scale_fill_manual(values = col_vals)+
  stat_summary(fun = median, geom = "point", colour = "white", size = 1) +
  #facet_wrap(~ thres, scales = "free_y") +
  
  facet_grid(fraction ~ thres)+
  labs(
    x = "Number of reads pairs",
    y = "Peaks per gene (log10 scale)",
    title = "B"
  ) +
  pub_theme +
  theme(axis.text = element_text(size = 10),
    legend.position = "none"
  )+ scale_x_discrete(
    labels = function(x) paste0(x, "M")
  )

legend_shared <- cowplot::get_legend(
  fig5a +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
)

# Remove legends from individual panels

fig5a_clean <- fig5a + theme(legend.position = "none")
fig5b_clean <- fig5b2 + theme(legend.position = "none")


# Combine top row

combined_fig5 <- 
  fig5a_clean + fig5b_clean + legend_shared +
  plot_layout(
    heights = c(1, 1.5, 0.15), nrow = 3
  )
ggsave(
  filename = "SupplementaryFig_5_dataset1.pdf",
  plot = combined_fig5,
  width = 12,
  height = 12
)

ggsave(
  filename = "SupplementaryFig_5_dataset1.png",
  plot = combined_fig5,
  width = 12,
  height = 12
)
