
##########################
##       Figure 6       ##
##########################

library(pheatmap)
library(RColorBrewer)
library(tibble)
library(dplyr)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)
library(tidyr)

###############
## functions ##
###############

genome_size <- 2.7e9   # hg38 mappable genome (approx)
## calculate forbes coefficient between two peak sets
forbes_coefficient <- function(gr1, gr2, genome_size ) {
  
  gr1 <- GenomicRanges::reduce(gr1)
  gr2 <- GenomicRanges::reduce(gr2)
  
  # force numeric to avoid overflow
  len1 <- sum(as.numeric(width(gr1)))
  len2 <- sum(as.numeric(width(gr2)))
  
  overlap <- GenomicRanges::intersect(gr1, gr2)
  len_overlap <- sum(as.numeric(width(overlap)))
  
  genome_size <- as.numeric(genome_size)
  
  (len_overlap * genome_size) / (len1 * len2)
}
## make heatmap plot of similarity between peak sets
make_plot_sim <- function(res, value_col, title, option_c,title.text,sub) {
  p <- ggplot(
    res,
    aes(
      x = group1_dataset1,
      y = group2_dataset2,
      fill = .data[[value_col]]
    )
  ) +
    geom_tile(
      colour = "grey95",
      linewidth = 0.25
    ) +
    geom_text(
      aes(
        label = sprintf("%.2f", .data[[value_col]])
      ),
      size = 6,
      colour = "black"
    ) +
    scale_fill_viridis_c(
      option = option_c,
      name = title,
      begin = 0.3,
      end = 0.9,
      na.value = "grey90",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(7, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +  
    coord_equal(expand = FALSE) +
    labs(
      x = "Dataset 1",
      y = "Dataset 2", title = title.text, subtitle = sub
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      plot.subtitle = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6), 
      plot.title = element_text(size = 24, face = "bold",hjust = -0.25)
    )
  return(p)
}
## simplify genomic annotation
changeAnnotation <- function(df_sub){
  annotation_df <- df_sub %>%
    dplyr::mutate(
      annotation_group = dplyr::case_when(
        grepl("Promoter", annotation) ~ "Promoter",
        grepl("Exon", annotation) ~ "Exon",
        grepl("Intron", annotation) ~ "Intron",
        TRUE ~ "Intergenic"
      )
    ) 
}
## simplify genomic annotation and output % of peaks in each category
changeAnnotation_perPeaks <- function(df_sub){
  annotation_df <- df_sub %>%
    dplyr::mutate(
      annotation_group = dplyr::case_when(
        grepl("Promoter", annotation) ~ "Promoter",
        grepl("Exon", annotation) ~ "Exon",
        grepl("Intron", annotation) ~ "Intron",
        TRUE ~ "Intergenic"
      )
    ) %>%
    dplyr::group_by(annotation_group) %>%
    dplyr::summarise(
      annotation_n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      percent = annotation_n / sum(annotation_n) * 100
    )
}

## overlap peaks (shared)
findPeaksOver <- function(gr1, gr2){
  
  hits <- findOverlaps(gr1, gr2, ignore.curatedand = TRUE)
  
  unique_query <- unique(queryHits(hits))
  
  return(gr1[unique_query])
}
## non-overlap peaks (specific)
findPeaksOverNo <- function(gr1,gr2){
  
  hits <- subsetByOverlaps(x = gr1, ranges = gr2, invert = TRUE)
  return(hits)
}
## per-peak-pair jaccard index
perPeakPairJaccard <- function(gr1, gr2) {
  
  hits <- findOverlaps(gr1, gr2, ignore.curatedand = TRUE)
  q <- queryHits(hits)
  s <- subjectHits(hits)
  overlap_bp <- width(pintersect(gr1[q], gr2[s]))
  peak1_width <- width(gr1[q])
  peak2_width <- width(gr2[s])
  
  pct_peak1 <- 100 * overlap_bp / peak1_width
  pct_peak2 <- 100 * overlap_bp / peak2_width
  
  jaccard <- overlap_bp /
    (peak1_width + peak2_width - overlap_bp)
  data.frame(
    peak1 = q,
    peak2 = s,
    peak1_width = peak1_width,
    peak2_width = peak2_width,
    overlap_bp = overlap_bp,
    pct_peak1 = pct_peak1,
    pct_peak2 = pct_peak2,
    jaccard = jaccard,
    jaccard_percent = 100 * jaccard
  )
}

## universal theme for plotting ##
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    curatedip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )

peakSets_fill_cols <- c(
  "NeuN default" = alpha("hotpink", 0.45),
  "NeuN curated" = "hotpink",
  "SOX10 default" = alpha("#6B8E23", 0.45),
  "SOX10 curated" = "#6B8E23")

######################
## Gather peak sets ##
######################

cf <- rep("NeuN+", "SOX10+",2)
curatedPeaks <-"./methodsPaper/data/curatedPeakSets/"
defaultPeaks <-"./methodsPaper/data/defaultedPeakSets/"

peakFiles<-c(paste0(defaultPeaks,"neun_q5e-2.narrowPeak.filt"),
                    paste0(defaultPeaks,"sox10_q5e-2.narrowPeak.filt"),
             paste0(curatedPeaks,"neun_q1e-16.narrowPeak.filt"),
             paste0(curatedPeaks,"sox10_q1e-18.narrowPeak.filt"))

print("Peaks files: ")
print(peakFiles)
peaks<-list()
for(i in 1:length(peakFiles)){
  peaks[[i]]<-makeGRangesFromDataFrame(read.table(file.path(peakFiles[i])), keep.extra.columns=TRUE, seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       ignore.curatedand=TRUE)
  peaks[[i]] <- peaks[[i]][!seqnames(peaks[[i]]) %in% c("chrX", "chrY"),]
  peaks[[i]] <- as.data.frame(peaks[[i]])
  colnames(peaks[[i]]) <- c("chr", "start", "end", "width", "curatedand", "peakID", "score", "curatedand2", "signal", "pvalue", "qvalue", "score2")
  
  peaks[[i]]$pvals <- 10^(-peaks[[i]]$pvalue)
  peaks[[i]]$qvals<- 10^(-peaks[[i]]$qvalue)
  peaks[[i]]$cf <- cf[i]
  peaks[[i]]$dataset <- "D2"
}
print(lapply(peaks, function(i) nrow(i)))
names(peaks) <- c("NeuN+_D2_default", "SOX10+_D2_default","NeuN+_D2_curated", "SOX10+_D2_curated")
peaks.D2.df <- peaks

peakFiles<-c(paste0(defaultPeaks,"NEUN_5e-2.narrowPeak.filt"),
             paste0(defaultPeaks,"SOX_5e-2.narrowPeak.filt"),
             paste0(curatedPeaks,"NEUN_1e-32.narrowPeak.filt"),
             paste0(curatedPeaks,"SOX_1e-25.narrowPeak.filt"))
print("Peaks files: ")
print(peakFiles)
peaks<-list()
for(i in 1:length(peakFiles)){
  peaks[[i]]<-makeGRangesFromDataFrame(read.table(file.path(paste0(peakFiles[i]))), keep.extra.columns=TRUE, seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       ignore.curatedand=TRUE)
  peaks[[i]] <- peaks[[i]][!seqnames(peaks[[i]]) %in% c("chrX", "chrY"),]
  peaks[[i]] <- as.data.frame(peaks[[i]])
  colnames(peaks[[i]]) <- c("chr", "start", "end", "width", "curatedand", "peakID", "score", "curatedand2", "signal", "pvalue", "qvalue", "score2")
  
  peaks[[i]]$pvals <- 10^(-peaks[[i]]$pvalue)
  peaks[[i]]$qvals<- 10^(-peaks[[i]]$qvalue)
  peaks[[i]]$cf <- cf[i]
  peaks[[i]]$dataset <- "D1"
}

print(lapply(peaks, function(i) nrow(i)))
names(peaks) <- c("NeuN+_D1_default", "SOX10+_D1_default","NeuN+_D1_curated", "SOX10+_D1_curated")
peaks.D1.df <- peaks

peak_groups.D1.gr <- lapply(peaks.D1.df, GRanges)
names(peak_groups.D1.gr) <- names(peaks.D1.df)
peak_groups.D2.gr <- lapply(peaks.D2.df, GRanges)
names(peak_groups.D2.gr) <- names(peaks.D2.df)
setPeaks <- c(names(peaks.D1.df), names(peaks.D2.df) )

## calculate forbes coefficient across peak sets ##
forbes_matrix <- matrix(data = NA, nrow =4, ncol =4)
for(i in 1:length(peak_groups.D1.gr)){
  for(j in 1:length(peak_groups.D2.gr)){
    forbes_matrix[i,j] <- forbes_coefficient(peak_groups.D1.gr[[i]], peak_groups.D2.gr[[j]],genome_size)
    forbes_matrix[j,i] <- forbes_coefficient(peak_groups.D2.gr[[j]], peak_groups.D1.gr[[i]],genome_size)
    
  }
}
colnames(forbes_matrix)<-names(peaks.D2.df)
rownames(forbes_matrix)<-names(peaks.D1.df)

forbes.plot <- as.data.frame(forbes_matrix) %>%
  rownames_to_column("cell1") %>%
  pivot_longer(
    -cell1,
    names_to = "cell2",
    values_to = "forbes"
  )

forbes.plot <- forbes.plot %>%
  mutate(
    group1_dataset1 = gsub("_(D1|D2)_", " ", cell1),
    group2_dataset2 = gsub("_(D1|D2)_", " ", cell2)
  )
forbes.plot$group1_dataset1 <- sub( "_", " ",forbes.plot$group1_dataset1) 

# ########################## Fig 6.a1 ############################
# Heatmap of forbes coefficient across all peaks in peak sets

p6a1 <- make_plot_sim(
  forbes.plot,
  "forbes",
  "Forbes coefficient",
  "plasma", "A","All peaks"
)

# ########################## Fig 6.a2 ############################
# Heatmap of forbes coefficient across only peaks annotated as 
# promoters in peak sets
             
allPeaks <- c(peak_groups.D1.gr, peak_groups.D2.gr)

peaksAnno.all<-lapply(allPeaks,annotatePeak,tssRegion=c(-3000, 3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", verbose = FALSE)

peaksAnno.df <- lapply(peaksAnno.all, function(x) as.data.frame(x@anno))
peaksAnno.df.anno <- lapply(peaksAnno.df, changeAnnotation)
peaksAnno.df.prom <- lapply(peaksAnno.df.anno, function(x) x[x$annotation_group == "Promoter",])

peaksAnno.df.prom.D1 <- peaksAnno.df.prom[c(1:4)]
peaksAnno.df.prom.D1 <- lapply(peaksAnno.df.prom.D1, GRanges)
names(peaksAnno.df.prom.D1) <- names(peaksAnno.df.prom.D1)

peaksAnno.df.prom.D2 <- peaksAnno.df.prom[c(5:8)]
peaksAnno.df.prom.D2 <- lapply(peaksAnno.df.prom.D2, GRanges)
names(peaksAnno.df.prom.D2) <- names(peaksAnno.df.prom.D2)

forbes_matrix.prom <- matrix(data = NA, nrow =4, ncol =4)
for(i in 1:length(peaksAnno.df.prom.D1)){
  for(j in 1:length(peaksAnno.df.prom.D2)){
    forbes_matrix.prom[i,j] <- forbes_coefficient(peaksAnno.df.prom.D1[[i]], peaksAnno.df.prom.D2[[j]],genome_size)
    forbes_matrix.prom[j,i] <- forbes_coefficient(peaksAnno.df.prom.D2[[j]], peaksAnno.df.prom.D1[[i]],genome_size)
    
  }
}
colnames(forbes_matrix.prom)<-names(peaksAnno.df.prom.D2)
rownames(forbes_matrix.prom)<-names(peaksAnno.df.prom.D1)

forbes.plot.prom <- as.data.frame(forbes_matrix.prom) %>%
  rownames_to_column("cell1") %>%
  pivot_longer(
    -cell1,
    names_to = "cell2",
    values_to = "forbes"
  )
forbes.plot.prom <- forbes.plot.prom %>%
  mutate(
    group1_dataset1 = gsub("_(D1|D2)_", " ", cell1),
    group2_dataset2 = gsub("_(D1|D2)_", " ", cell2)
  )

p6a2 <- make_plot_sim(
  forbes.plot.prom,
  "forbes",
  "Forbes coefficient",
  "plasma", "", "Promoter peaks"
)

## Peaks that overlap vs no overlap
NEUNover_default.D1 <- findPeaksOver(peak_groups.D1.gr[[1]], peak_groups.D2.gr[[1]])
NEUNover_default.D2 <- findPeaksOver(peak_groups.D2.gr[[1]], peak_groups.D1.gr[[1]])

NEUNoverno_default.D1 <- findPeaksOverNo(peak_groups.D1.gr[[1]], peak_groups.D2.gr[[1]])
NEUNoverno_default.D2 <- findPeaksOverNo(peak_groups.D2.gr[[1]], peak_groups.D1.gr[[1]])

NEUNover_curated.D1 <- findPeaksOver(peak_groups.D1.gr[[3]], peak_groups.D2.gr[[3]])
NEUNover_curated.D2 <- findPeaksOver(peak_groups.D2.gr[[3]], peak_groups.D1.gr[[3]])

NEUNoverno_curated.D1 <- findPeaksOverNo(peak_groups.D1.gr[[3]], peak_groups.D2.gr[[3]])
NEUNoverno_curated.D2 <- findPeaksOverNo(peak_groups.D2.gr[[3]], peak_groups.D1.gr[[3]])

overPeaks.NEUN <- list(NEUNover_default.D1,NEUNover_default.D2,
                  NEUNoverno_default.D1,NEUNoverno_default.D2,
                  NEUNover_curated.D1,NEUNover_curated.D2,
                  NEUNoverno_curated.D1,NEUNoverno_curated.D2)
names(overPeaks.NEUN)<- c("NeuN+\nshared\nQ5e-2\nD1", "NeuN+\nnshared\nQ5e-2\nD2",
                     "NeuN+\nspecific\nQ5e-2\nD1", "NeuN+\nspecific\nQ5e-2\nD2",
                     "NeuN+\nnshared\nQ1e-32\nD1", "NeuN+\nnshared\nQ1e-16\nD2",
                     "NeuN+\nspecific\nQ1e-32\nD1", "NeuN+\nspecific\nQ1e-16\nD2")
peaksAnno.all.NEUN<-lapply(overPeaks.NEUN,annotatePeak,
                      TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                      annoDb="org.Hs.eg.db", verbose = FALSE)
names(peaksAnno.all.NEUN) <- names(overPeaks.NEUN)

peaksAnno.df.NEUN <- lapply(peaksAnno.all.NEUN, function(x) as.data.frame(x@anno))
peaksAnno.df.anno.NEUN <- lapply(peaksAnno.df.NEUN, changeAnnotation_perPeaks)

dataset <- rep(c("Dataset 1","Dataset 2"), 4)
mode <- c(rep(c("shared", "specific"),each =2),rep(c("shared", "specific"),each =2))
thres <- c(rep("default", 4),rep("curated",4))
peaksets <- rep(c("Dataset 1","Dataset 2"), 4)
peaksAnno.all.df.NEUN <- lapply(seq_along(peaksAnno.df.anno.NEUN), function(i) {
  df <- peaksAnno.df.anno.NEUN[[i]]
  df$mode <- mode[i]
  df$dataset <- dataset[i]
  df$thres <- thres[i]
  df$peakset <- paste0(peaksets[i]," ",mode[i]," ",thres[i])
  df
})
peaksAnno.all.df.NEUN <- bind_rows(peaksAnno.all.df.NEUN)
peaksAnno.all.df.NEUN$peakset <- as.character(peaksAnno.all.df.NEUN$peakset)
level_order <- c("Dataset 1 shared default", "Dataset 1 specific default",
                 "Dataset 2 shared default", "Dataset 2 specific default",
                 "Dataset 1 shared curated", "Dataset 1 specific curated",
                 "Dataset 2 shared curated", "Dataset 2 specific curated" )
peaksAnno.all.df.NEUN$peakset <- factor(peaksAnno.all.df.NEUN$peakset,
                                   levels =rev(level_order ))
                            
# ########################## Fig 6.b1 ############################
# Stacked bar plot of % of peaks per genomic annotation (x-axis)
# for each peak set (NEUN only)
                            
p6c1 <- ggplot(
  peaksAnno.all.df.NEUN,
  aes(
    x = peakset,
    y = percent,
    fill = annotation_group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  labs(
    x = "",
    y = "% Peaks",
    fill = "Annotation", title = "C", subtitle = "NeuN+"
  ) +coord_flip()+ 
  pub_theme +theme(
    plot.title = element_text(size = 24,face = "bold",hjust = -0.5),
    plot.subtitle = element_text(size = 22,
                                 face = "bold",
                                 hjust = 0.5
    ))+ guides(fill = guide_legend(ncol = 2))
                            

## Peaks that overlap vs no overlap
SOXover_default.D1 <- findPeaksOver(peak_groups.D1.gr[[2]], peak_groups.D2.gr[[2]])
SOXover_default.D2 <- findPeaksOver(peak_groups.D2.gr[[2]], peak_groups.D1.gr[[2]])

SOXoverno_default.D1 <- findPeaksOverNo(peak_groups.D1.gr[[2]], peak_groups.D2.gr[[2]])
SOXoverno_default.D2 <- findPeaksOverNo(peak_groups.D2.gr[[2]], peak_groups.D1.gr[[2]])

SOXover_curated.D1 <- findPeaksOver(peak_groups.D1.gr[[4]], peak_groups.D2.gr[[4]])
SOXover_curated.D2 <- findPeaksOver(peak_groups.D2.gr[[4]], peak_groups.D1.gr[[4]])

SOXoverno_curated.D1 <- findPeaksOverNo(peak_groups.D1.gr[[4]], peak_groups.D2.gr[[4]])
SOXoverno_curated.D2 <- findPeaksOverNo(peak_groups.D2.gr[[4]], peak_groups.D1.gr[[4]])

overPeaksSOX <- list(SOXover_default.D1,SOXover_default.D2,
                  SOXoverno_default.D1,SOXoverno_default.D2,
                  SOXover_curated.D1,SOXover_curated.D2,
                  SOXoverno_curated.D1,SOXoverno_curated.D2)
peaksAnno.all.SOX<-lapply(overPeaksSOX,annotatePeak,
                      TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                      annoDb="org.Hs.eg.db", verbose = FALSE)
names(peaksAnno.all.SOX) <- names(overPeaksSOX)
peaksAnno.df.SOX <- lapply(peaksAnno.all.SOX, function(x) as.data.frame(x@anno))
peaksAnno.df.anno.SOX <- lapply(peaksAnno.df.SOX, changeAnnotation)
peaksAnno.all.df.SOX <- lapply(seq_along(peaksAnno.df.anno.SOX), function(i) {
  df <- peaksAnno.df.anno.SOX[[i]]
  df$mode <- mode[i]
  df$dataset <- dataset[i]
  df$thres <- thres[i]
  df$peakset <- paste0(peaksets[i]," ",mode[i]," ",thres[i])
  df
})
peaksAnno.all.df.SOX <- bind_rows(peaksAnno.all.df.SOX)
peaksAnno.all.df.SOX$peakset <- as.character(peaksAnno.all.df.SOX$peakset)
level_order <- c("Dataset 1 shared default", "Dataset 1 specific default",
                 "Dataset 2 shared default", "Dataset 2 specific default",
                 "Dataset 1 shared curated", "Dataset 1 specific curated",
                 "Dataset 2 shared curated", "Dataset 2 specific curated" )
peaksAnno.all.df.SOX$peakset <- factor(peaksAnno.all.df.SOX$peakset,
                                   levels =rev(level_order ))
                           
# ########################## Fig 6.b2 ############################
# Stacked bar plot of % of peaks per genomic annotation (x-axis)
# for each peak set (SOX10 only)

p6c2<- ggplot(
  peaksAnno.all.df.SOX,
  aes(
    x = peakset,
    y = percent,
    fill = annotation_group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  labs(
    x = "",
    y = "% Peaks",
    fill = "Annotation", subtitle = "SOX10+"
  ) +coord_flip()+ 
  pub_theme +theme(
    legend.position = "right",
    plot.title = element_text(size = 24,face = "bold"),
    plot.subtitle = element_text(size = 22,
                                 face = "bold",
                                 hjust = 0.5
    ))+ guides(fill = guide_legend(ncol = 1))

                           
## calculate per-peak-pair Jaccard scores
D1_D2_default_NEUN <- perPeakPairJaccard(
  peak_groups.D1.gr[[1]],
  peak_groups.D2.gr[[1]]
)

D1_D2_curated_NEUN <- perPeakPairJaccard(
  peak_groups.D1.gr[[3]],
  peak_groups.D2.gr[[3]]
)
D1_D2_default_SOX <- perPeakPairJaccard(
  peak_groups.D1.gr[[2]],
  peak_groups.D2.gr[[2]]
)

D1_D2_curated_SOX <- perPeakPairJaccard(
  peak_groups.D1.gr[[4]],
  peak_groups.D2.gr[[4]]
)

D1_D2_default_NEUN$threshold <- "default"
D1_D2_curated_NEUN$threshold <- "curated"
D1_D2_default_NEUN$cellType <- "NeuN+"
D1_D2_curated_NEUN$cellType <- "NeuN+"
D1_D2_default_SOX$threshold <- "default"
D1_D2_curated_SOX$threshold <- "curated"
D1_D2_default_SOX$cellType <- "SOX10+"
D1_D2_curated_SOX$cellType <- "SOX10+"
D1_D2$peaksets <- paste0(D1_D2$cellType, " ",D1_D2$threshold )

# ########################## Fig 6.b ############################
# Violin plot of per-peak-pair Jaccard scores (y-axis)
# per peak set comparison
# coloured by cell type (NEUN/SOX10) and peak set (curated/default)
                           
p6b <- ggplot(
  D1_D2,
  aes(
    x = peaksets,
    y = jaccard_percent, fill = peaksets )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.10,
    fill = "black",
    outlier.shape = NA
  ) + stat_summary(
    fun = median,position = position_dodge(width = 0.7),
    geom = "point",
    size = 2.5,
    color = "white"
  )+
  scale_fill_manual(values = peakSets_fill_cols)+
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    x = "Comparison of peak sets",
    y = "Per-peak-pair Jaccard overlap (%)", title = "B"
  ) +
  pub_theme +  theme_bw(12)+
  theme(axis.title = element_text(size = 14), legend.position = "none",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 24,face = "bold"),
        plot.subtitle = element_text(size = 22,
                                     face = "bold",
                                     hjust = 0.5))

## gather plots for Figure 6 ##
panel_top <- (p6a1 | p6a2) +
  plot_layout(widths = c(1.5,1.5))

p6c1_clean <- p6c1 + theme(legend.position = "none")
p6c2_clean <- p6c2 + theme(legend.position = "right")

bottom_pannel <- p6b | p6c1_clean | p6c2_clean |  
                           
figure6 <- panel_top / bottom_pannel  +
  plot_layout(heights = c(1.25,1))

ggsave(
  filename = "Figure6.png",
  plot = figure6,
  width = 20,
  height = 15
)

ggsave(
  filename = "Figure6.pdf",
  plot = figure6,
  width = 20,
  height = 15
)
