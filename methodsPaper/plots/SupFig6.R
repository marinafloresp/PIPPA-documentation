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


groupsPeaks <-"/path/to/peaks/"
brainGenes <- read.table("brainGenes.txt")[,1]
brainGenes <- unique(brainGenes)
brainGenes <- sub("\\..*$", "", brainGenes)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")
peakFiles<-list.files(groupsPeaks, pattern="_1e-2.*\\.narrowPeak\\.filt|_5e-2.*\\.narrowPeak\\.filt|_1e-5.*\\.narrowPeak\\.filt")
files_readsAll <- peakFiles[grepl("_\\d+readsAll\\.narrowPeak\\.filt", peakFiles)]
files_readsAll <- files_readsAll[grepl("SOX|NEUN", files_readsAll, ignore.case = TRUE)]
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

peak_sets_sub <- peaks[grepl("5e-2", names(peaks))]
peakAnno.all <- lapply(peaks, annotatePeak,
                       TxDb = txdb,
                       annoDb = "org.Hs.eg.db")
names(peakAnno.all) <- names(peaks)
peakAnno.all <- peakAnno.all[grepl("SOX|NEUN", names(peakAnno.all), ignore.case = TRUE)]
peakAnno.all <- peakAnno.all[grepl("_\\d+readsAll", names(peakAnno.all))]
calculatePeaksGenes <- function(peakAnno, setName){
  annot <- peakAnno %>%
    as.data.frame()
  
  peaks_per_gene <- annot %>%
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
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "SOX"] <- "SOX10+"
genesReadsSets.all$thres <-fct_recode(genesReadsSets.all$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")


# ########################## Fig Sup 6.a ###########################################
# Line graph of number of reads included (x-axis)
# by number of total genes annotated (y-axis)
# coloured by q-value threshold (color) 
# line-type by cell-type

fig6a <- ggplot(as.data.frame(genesReadsSets.all), aes(x = num_reads, y = n_genes, color = thres, linetype = fraction)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSets.all$num_reads))
  ) +
  scale_linetype_manual(values = c("NeuN+" = "solid", "SOX10+" = "dotted")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Read pairs (Millions)",
    y = "Total annotated genes",
    title = "A",
    color = "Q value threshold", linetype = "Cell type"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 6.b ###########################################
# Line graph of number of reads included (x-axis)
# by mean number of peaks annotated to genes (y-axis)
# coloured by q-value threshold (color) 
# line-type by cell-type

fig6b <- ggplot(
  genesReadsSets.all,
  aes(x = factor(num_reads), y = mean_peaks_per_gene, color = thres,
      linetype = fraction,
      group = interaction(thres, fraction))
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_discrete(
    breaks = sort(unique(genesReadsSets.all$num_reads))
  )+
  # scale_y_log10()+
  scale_linetype_manual(values = c("NeuN+" = "solid", "SOX10+" = "dotted")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Read pairs (Millions)",
    y = "Mean peaks per gene", 
    color = "Q value threshold", 
    linetype = "Cell type",
    title = "B"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )


# ########################## Fig Sup 6.c ###########################################
# Line graph of number of reads included (x-axis)
# by % of peaks annotated as brain-related genes
# coloured by q-value threshold (color) 
# line-type by cell-type

calculateBrainGenes <- function(peakAnno, setName){
  annot <- peakAnno@anno %>%
    as.data.frame()
  
  annot$ENSEMBL <- sub("\\..*$", "", annot$ENSEMBL)
  
  ids <- unique(annot$ENSEMBL)
  genes_with_peaks <- ids
  intersect <- intersect(genes_with_peaks,brainGenes)
  
  perBrainGenes <- length(intersect) / length(genes_with_peaks) 
  
  list(
    summary = tibble::tibble(
      name = setName,
      propBrainGenes = perBrainGenes,
      
    )
  )
}

genesReadsSetsBrain <- Map(
  calculateBrainGenes,
  peakAnno.all,
  names(peakAnno.all)
)

all_summaries.Brain <- lapply(genesReadsSetsBrain, `[[`, "summary")
genesReadsSetsBrain.all <- bind_rows(all_summaries.Brain)
genesReadsSetsBrain.all$num_reads <- as.integer(str_extract(genesReadsSetsBrain.all$name, "(?<=_)[0-9]+(?=readsAll)"))
genesReadsSetsBrain.all$thres <- str_extract(genesReadsSetsBrain.all$name, "[0-9]+e-[0-9]+")
genesReadsSetsBrain.all$thres <-fct_recode(genesReadsSetsBrain.all$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")

genesReadsSetsBrain.all$thres <- factor(
  genesReadsSetsBrain.all$thres,
  levels = c("0.05" ,"0.01" , "0.00001")
)
genesReadsSetsBrain.all$fraction <-str_extract(genesReadsSetsBrain.all$name, "^[^_]+")
genesReadsSetsBrain.all$fraction[genesReadsSetsBrain.all$fraction == "NEUN"] <- "NeuN+"
genesReadsSetsBrain.all$fraction[genesReadsSetsBrain.all$fraction == "SOX"] <- "SOX10+"

fig6c <- ggplot(
  genesReadsSetsBrain.all,
  aes(
    x = num_reads,
    y = propBrainGenes*100,
    color = thres,
    linetype = fraction,
    group = interaction(thres, fraction)
  )
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSetsBrain.all$num_reads))
  ) +
  scale_linetype_manual(values = c("NeuN+" = "solid", "SOX10+" = "dotted")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Read pairs (Millions)",
    y = "% Peaks annotated as brain genes",
    color = "Q value threshold", linetype = "Cell type",
    title = "C"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 6.d ###########################################
# Line graph of number of reads included (x-axis)
# by mean number of narrow peaks annotated to genes
# coloured by q-value threshold (color) 
# line-type by cell-type

calculatePeaksGenes <- function(peakAnno, setName){
  
  annot <- peakAnno %>%
    as.data.frame() %>%
    dplyr::mutate(
      peak_size = end - start + 1,
      short_peak = peak_size < 400,
      promoter_peak = grepl("Promoter", annotation)
    )
  
  peaks_per_gene <- annot %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      
      # All peaks
      peak_count = dplyr::n(),
      
      # Peaks <400 bp
      peaks_lt400 = sum(short_peak),
      
      # Promoter peaks, regardless of size
      promoter_peaks = sum(promoter_peak),
      
      # Fraction of peaks <400 bp
      fraction_lt400 = peaks_lt400 / peak_count,
      
      # Fraction of peaks that are promoters
      fraction_promoter = promoter_peaks / peak_count,
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      peak_set = setName
    )
  
  list(
    summary = tibble::tibble(
      
      name = setName,
      
      # Genes
      n_genes = nrow(peaks_per_gene),
      
      # All peaks
      mean_peaks_per_gene = mean(peaks_per_gene$peak_count),
      total_annotated_peaks = sum(peaks_per_gene$peak_count),
      
      # <400 bp peaks
      total_peaks_lt400 = sum(peaks_per_gene$peaks_lt400),
      fraction_peaks_lt400 =
        sum(peaks_per_gene$peaks_lt400) /
        sum(peaks_per_gene$peak_count),
      mean_peaks_lt400 =
        mean(peaks_per_gene$peaks_lt400),
      
      # Promoter peaks
      total_promoter_peaks =
        sum(peaks_per_gene$promoter_peaks),
      mean_promoter_peaks =
        mean(peaks_per_gene$promoter_peaks),
      fraction_promoter_peaks =
        sum(peaks_per_gene$promoter_peaks) /
        sum(peaks_per_gene$peak_count)
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
genesReadsSetsNarrow <- bind_rows(all_summaries)
genesReadsSetsNarrow$num_reads <- as.integer(str_extract(genesReadsSetsNarrow$name, "(?<=_)[0-9]+(?=readsAll)"))
genesReadsSetsNarrow$thres <- str_extract(genesReadsSetsNarrow$name, "[0-9]+e-[0-9]+")
genesReadsSetsNarrow$thres <-fct_recode(genesReadsSetsNarrow$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")

genesReadsSetsNarrow$thres <- factor(
  genesReadsSetsNarrow$thres,
  levels = c("0.05" ,"0.01" , "0.00001")
)
genesReadsSetsNarrow$fraction <-str_extract(genesReadsSetsNarrow$name, "^[^_]+")
genesReadsSetsNarrow$fraction[genesReadsSetsNarrow$fraction == "NEUN"] <- "NeuN+"
genesReadsSetsNarrow$fraction[genesReadsSetsNarrow$fraction == "SOX"] <- "SOX10+"
genesReadsSetsNarrow[genesReadsSetsNarrow$num_reads == 50 & genesReadsSetsNarrow$thres == "0.05",]$mean_peaks_lt400

fig6d <- ggplot(
  genesReadsSetsNarrow,
  aes(
    x = num_reads,
    y =mean_peaks_lt400,
    color = thres,
    linetype = fraction,
    group = interaction(thres, fraction)
  )
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSetsNarrow$num_reads))
  ) +
  scale_linetype_manual(values = c("NeuN+" = "solid", "SOX10+" = "dotted")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Read pairs (Millions)",
    y = "Mean narrow peaks (<400bp) per gene",
    color = "Q value threshold", linetype = "Cell type",
    title = "D"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 6.e ###########################################
# Line graph of number of reads included (x-axis)
# by mean number of peaks annotated as promoters per gene
# coloured by q-value threshold (color) 
# line-type by cell-type

fig6e <- ggplot(
  genesReadsSetsNarrow,
  aes(
    x = num_reads,
    y =mean_promoter_peaks,
    color = thres,
    linetype = fraction,
    group = interaction(thres, fraction)
  )
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSetsNarrow$num_reads))
  ) +
  scale_linetype_manual(values = c("NeuN+" = "solid", "SOX10+" = "dotted")) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = "Read pairs (Millions)",
    y = "Mean promoter peaks per gene",
    color = "Q value threshold", linetype = "Cell type",
    title = "E"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

## Gather plots for Supplementary Figure 6 ##


fig6a_clean <- fig6a + theme(legend.position = "none")
fig6b_clean <- fig6b + theme(legend.position = "none")
fig6c_clean <- fig6c + theme(legend.position = "none")
fig6d_clean <- fig6d + theme(legend.position = "none")
fig6e_clean <- fig6e + theme(legend.position = "none")
legend_shared <- cowplot::get_legend(
  fig6a +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
)

# Combine top row

combined_Supfig6 <- 
  (fig6a_clean | fig6b_clean) / (fig6c_clean | fig6d_clean) / (fig6e_clean | legend_shared) +
  plot_layout(
    heights = c(1, 1, 1), nrow = 3
  )
combined_Supfig6

ggsave(
  filename = "SupplementaryFig_6_dataset1.pdf",
  plot = combined_Supfig6,
  width = 14,
  height = 16
)

ggsave(
  filename = "SupplementaryFig_6_dataset1.png",
  plot = combined_Supfig6,
  width = 14,
  height = 16
)


