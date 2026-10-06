#####################################
##       Supplementary Fig 14      ##
#####################################

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

## overlap vs non-overlap peaks 
findPeaksOver <- function(gr1, gr2){
  
  hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  
  unique_query <- unique(queryHits(hits))
  
  return(gr1[unique_query])
}

findPeaksOverNo <- function(gr1,gr2){
  
  hits <- subsetByOverlaps(x = gr1, ranges = gr2, invert = TRUE)
  return(hits)
}

## theme for plotting
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )

colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")

######################
## Gather peak sets ##
######################

cf <- rep("NeuN+", "SOX10+",2)
strPeaks <-"./methodsPaper/data/strPeakSets/"
relaxPeaks <-"./methodsPaper/data/relaxedPeakSets/"

peakFiles<-c(paste0(relaxPeaks,"neun_q5e-2.narrowPeak.filt"),
             paste0(relaxPeaks,"sox10_q5e-2.narrowPeak.filt"),
             paste0(strPeaks,"neun_q1e-16.narrowPeak.filt"),
             paste0(strPeaks,"sox10_q1e-18.narrowPeak.filt"))

print("Peaks files: ")
print(peakFiles)
peaks<-list()
for(i in 1:length(peakFiles)){
  peaks[[i]]<-makeGRangesFromDataFrame(read.table(file.path(peakFiles[i])), keep.extra.columns=TRUE, seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       ignore.strand=TRUE)
  peaks[[i]] <- peaks[[i]][!seqnames(peaks[[i]]) %in% c("chrX", "chrY"),]
  peaks[[i]] <- as.data.frame(peaks[[i]])
  colnames(peaks[[i]]) <- c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal", "pvalue", "qvalue", "score2")
  
  peaks[[i]]$pvals <- 10^(-peaks[[i]]$pvalue)
  peaks[[i]]$qvals<- 10^(-peaks[[i]]$qvalue)
  peaks[[i]]$cf <- cf[i]
  peaks[[i]]$dataset <- "D2"
}
print(lapply(peaks, function(i) length(i)))
names(peaks) <- c("NeuN+_D2_default", "SOX10+_D2_default","NeuN+_D2_curated", "SOX10+_D2_curated")
peaks.D2.df <- peaks

peakFiles<-c(paste0(relaxPeaks,"NEUN_5e-2.narrowPeak.filt"),
             paste0(relaxPeaks,"SOX_5e-2.narrowPeak.filt"),
             paste0(strPeaks,"NEUN_1e-32.narrowPeak.filt"),
             paste0(strPeaks,"SOX_1e-25.narrowPeak.filt"))
print("Peaks files: ")
print(peakFiles)
peaks<-list()
for(i in 1:length(peakFiles)){
  peaks[[i]]<-makeGRangesFromDataFrame(read.table(file.path(paste0(peakFiles[i]))), keep.extra.columns=TRUE, seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       ignore.strand=TRUE)
  peaks[[i]] <- peaks[[i]][!seqnames(peaks[[i]]) %in% c("chrX", "chrY"),]
  peaks[[i]] <- as.data.frame(peaks[[i]])
  colnames(peaks[[i]]) <- c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal", "pvalue", "qvalue", "score2")
  
  peaks[[i]]$pvals <- 10^(-peaks[[i]]$pvalue)
  peaks[[i]]$qvals<- 10^(-peaks[[i]]$qvalue)
  peaks[[i]]$cf <- cf[i]
  peaks[[i]]$dataset <- "D1"
}

print(lapply(peaks, function(i) length(i)))
names(peaks) <- c("NeuN+_D1_default", "SOX10+_D1_default","NeuN+_D1_curated", "SOX10+_D1_curated")
peaks.D1.df <- peaks

peak_groups.D1.gr <- lapply(peaks.D1.df, GRanges)
names(peak_groups.D1.gr) <- names(peaks.D1.df)
peak_groups.D2.gr <- lapply(peaks.D2.df, GRanges)
names(peak_groups.D2.gr) <- names(peaks.D2.df)
setPeaks <- c(names(peaks.D1.df), names(peaks.D2.df) )

## Peaks that overlap vs no overlap for NEUN peak sets ##

NEUNover_relax.D1 <- findPeaksOver(peak_groups.D1.gr[[1]], peak_groups.D2.gr[[1]])
NEUNover_relax.D2 <- findPeaksOver(peak_groups.D2.gr[[1]], peak_groups.D1.gr[[1]])

NEUNoverno_relax.D1 <- findPeaksOverNo(peak_groups.D1.gr[[1]], peak_groups.D2.gr[[1]])
NEUNoverno_relax.D2 <- findPeaksOverNo(peak_groups.D2.gr[[1]], peak_groups.D1.gr[[1]])

NEUNover_str.D1 <- findPeaksOver(peak_groups.D1.gr[[3]], peak_groups.D2.gr[[3]])
NEUNover_str.D2 <- findPeaksOver(peak_groups.D2.gr[[3]], peak_groups.D1.gr[[3]])

NEUNoverno_str.D1 <- findPeaksOverNo(peak_groups.D1.gr[[3]], peak_groups.D2.gr[[3]])
NEUNoverno_str.D2 <- findPeaksOverNo(peak_groups.D2.gr[[3]], peak_groups.D1.gr[[3]])

overPeaks.NEUN <- list(NEUNover_relax.D1,NEUNover_relax.D2,
                  NEUNoverno_relax.D1,NEUNoverno_relax.D2,
                  NEUNover_str.D1,NEUNover_str.D2,
                  NEUNoverno_str.D1,NEUNoverno_str.D2)
names(overPeaks.NEUN)<- c("NeuN+\nshared\ndefault\nD1", "NeuN+\nnshared\ndefault\nD2",
                     "NeuN+\nspecific\ndefault\nD1", "NeuN+\nspecific\ndefault\nD2",
                     "NeuN+\nnshared\ncurated\nD1", "NeuN+\nnshared\ncurated\nD2",
                     "NeuN+\nspecific\ncurated\nD1", "NeuN+\nspecific\ncurated\nD2")
peaksAnno.all<-lapply(overPeaks.NEUN,annotatePeak,
                      TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                      annoDb="org.Hs.eg.db", verbose = FALSE)
names(peaksAnno.all) <- names(overPeaks.NEUN)
peaksAnno.df.NEUN <- lapply(peaksAnno.all, function(x) as.data.frame(x@anno))

## Peaks that overlap vs no overlap for SOX10 peak sets ##

SOXover_relax.D1 <- findPeaksOver(peak_groups.D1.gr[[2]], peak_groups.D2.gr[[2]])
SOXover_relax.D2 <- findPeaksOver(peak_groups.D2.gr[[2]], peak_groups.D1.gr[[2]])

SOXoverno_relax.D1 <- findPeaksOverNo(peak_groups.D1.gr[[2]], peak_groups.D2.gr[[2]])
SOXoverno_relax.D2 <- findPeaksOverNo(peak_groups.D2.gr[[2]], peak_groups.D1.gr[[2]])

SOXover_str.D1 <- findPeaksOver(peak_groups.D1.gr[[4]], peak_groups.D2.gr[[4]])
SOXover_str.D2 <- findPeaksOver(peak_groups.D2.gr[[4]], peak_groups.D1.gr[[4]])

SOXoverno_str.D1 <- findPeaksOverNo(peak_groups.D1.gr[[4]], peak_groups.D2.gr[[4]])
SOXoverno_str.D2 <- findPeaksOverNo(peak_groups.D2.gr[[4]], peak_groups.D1.gr[[4]])

overPeaksSOX <- list(SOXover_relax.D1,SOXover_relax.D2,
                     SOXoverno_relax.D1,SOXoverno_relax.D2,
                     SOXover_str.D1,SOXover_str.D2,
                     SOXoverno_str.D1,SOXoverno_str.D2)
names(overPeaksSOX)<- c("SOX10\nshared\ndefault\nD1", "SOX10\nnshared\ndefault\nD2",
                        "SOX10\nspecific\ndefault\nD1", "SOX10\nspecific\ndefault\nD2",
                        "SOX10\nnshared\ncurated\nD1", "SOX10\nnshared\ncurated\nD2",
                        "SOX10\nspecific\ncurated\nD1", "SOX10\nspecific\ncurated\nD2")
peaksAnno.all.SOX<-lapply(overPeaksSOX,annotatePeak,
                      TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                      annoDb="org.Hs.eg.db", verbose = FALSE)
names(peaksAnno.all.SOX) <- names(overPeaksSOX)
peaksAnno.df.SOX<- lapply(peaksAnno.all.SOX, function(x) as.data.frame(x@anno))

# ########################## Fig Sup 14.a ###########################################
# Stacked bar plots of % of peaks that are shared vs specific to curated and default peak sets

allPeaks <- c(overPeaks.NEUN,overPeaks.SOX)
nPeaks<-unlist(lapply(allPeaks, length))
totLengthPeaks<-unlist(lapply(lapply(allPeaks, width), sum))
avgLengthPeaks <- totLengthPeaks/nPeaks
allWidth <- unlist(lapply(allPeaks, width))

dataset <- rep(c("Dataset 1","Dataset 2"), 8)
mode <- c(rep(c("shared", "specific"),each =2),rep(c("shared", "specific"),each =2),
          rep(c("shared", "specific"),each =2),rep(c("shared", "specific"),each =2))
thres <- c(rep("default", 4),rep("curated",4),rep("default", 4),rep("curated",4))
peaksets <- rep(c("Dataset 1","Dataset 2"), 8)
ct <- rep(c("NeuN","SOX10"),each = 8)

numPeaks <- data.frame(
  peakset = paste0(peaksets," ",ct, " ",thres),
  over = mode, thres = thres,
  num = nPeaks,
  avglength = round(avgLengthPeaks,2),
  medlength = avgLengthPeaks, 
  dataset = dataset,
  ct = ct)
numPeaks <- numPeaks %>%
  group_by(dataset, peakset) %>%
  mutate(percent = num / sum(num) * 100) %>%
  ungroup()

numPeaks.NEUN <- numPeaks[numPeaks$ct == "NeuN",]

level_order <- c("Dataset 1 NeuN default", "Dataset 1 NeuN curated",
                 "Dataset 2 NeuN default", "Dataset 2 NeuN curated" )
numPeaks.NEUN$peakset <- factor(numPeaks.NEUN$peakset,
                                levels =rev(level_order ))
p14A1 <- ggplot(
  numPeaks.NEUN,
  aes(x = percent, y = peakset, fill = over)
) +
  geom_col(
    position = "stack",
    linewidth = 0.1,
    width = 0.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "shared" = "hotpink",
      "specific" = scales::alpha("hotpink", 0.45)
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "A", subtitle = "NeuN",
    y = "Peak set",
    x = "% of peaks",
    fill = " "
  ) +
  pub_theme +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14,
                                margin = margin(r = 0.1)),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 24, face = "bold", hjust = -0.20),
    plot.subtitle = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    )
  )

numPeaks.SOX <- numPeaks[numPeaks$ct == "SOX10",]
level_order <- c("Dataset 1 SOX10 default", "Dataset 1 SOX10 curated",
                 "Dataset 2 SOX10 default", "Dataset 2 SOX10 curated" )

numPeaks.SOX$peakset <- factor(numPeaks.SOX$peakset,
                               levels =rev(level_order ))

p14A2 <- ggplot(
  numPeaks.SOX,
  aes(x = percent, y = peakset, fill = over)
) +
  geom_col(
    position = "stack",
    linewidth = 0.1,
    width = 0.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "shared" = "#6B8E23",
      "specific" = scales::alpha("#6B8E23", 0.45)
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = " ", subtitle = "SOX10",
    y = "Peak set",
    x = "% of peaks",
    fill = " "
  ) +
  pub_theme +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14,
                                margin = margin(r = 0.1)),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 24, face = "bold", hjust = -0.15),
    plot.subtitle = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    )
  )
                          
# ########################## Fig Sup 14.b ###########################################
# Stacked bar plots of % of genes annotated to peaks that are related to brain for 
# shared vs specific peak sets                          

brainGenes <- read.table("./brainGenes.txt")[,1]
brainGenes <- unique(brainGenes)

brainGenes <- sub("\\..*$", "", brainGenes)
peaksAnno.df.all <- c(peaksAnno.df.NEUN, peaksAnno.df.SOX)
peaksAnno.df.all <- lapply(peaksAnno.df.all, function(df) {
  df$ENSEMBL <- sub("\\..*$", "", df$ENSEMBL)
  df
})

brainGenes.df <- lapply(seq_along(peaksAnno.df.all), function(i) {
  
  df <- peaksAnno.df.all[[i]]
  
  # Unique genes annotated to peaks
  genes_with_peaks <- unique(na.omit(df$ENSEMBL))
  
  # Brain-related genes
  brain_genes <- intersect(genes_with_peaks, brainGenes)
  
  total_genes <- length(genes_with_peaks)
  n_brain_genes <- length(brain_genes)
  
  data.frame(
    peakset = paste0(
      dataset[i], " ",
      mode[i], " ",
      thres[i]
    ),
    cellType = ct[i],
    total_genes = total_genes,
    brain_genes = n_brain_genes,
    non_brain_genes = total_genes - n_brain_genes
  )
})

brainGenes.df <- do.call(rbind, brainGenes.df)
rownames(brainGenes.df) <- NULL
brainGenes.df$peakset <- factor(brainGenes.df$peakset)

brainGenes.plot <- brainGenes.df %>%
  pivot_longer(
    cols = c(brain_genes, non_brain_genes),
    names_to = "annotation",
    values_to = "num"
  ) %>%
  mutate(
    annotation = recode(
      annotation,
      brain_genes = "Brain Gene",
      non_brain_genes = "Non-brain Gene"
    )
  )

brainGenes.plot$annotation <- factor(
  brainGenes.plot$annotation,
  levels = rev(c("Brain Gene", "Non-brain Gene"))
)

brainGenes.plot.NEUN <- brainGenes.plot[brainGenes.plot$cellType == "NeuN",]
level_order <- c("Dataset 1 specific default", "Dataset 1 shared default",
                 "Dataset 1 specific curated","Dataset 1 shared curated",
                 "Dataset 2 specific default", "Dataset 2 shared default",
                 "Dataset 2 specific curated","Dataset 2 shared curated")
brainGenes.plot.NEUN$peakset <- factor(brainGenes.plot.NEUN$peakset,
                                       levels =rev(level_order ))

p14B1 <- ggplot(
  brainGenes.plot.NEUN,
  aes(
    y = peakset,
    x = num,
    fill = annotation
  )
) +
  geom_col(
    position = "fill",
    width = 0.8,
    linewidth = 0,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) x * 100
  )+ 
  scale_fill_manual(
    values = c(
      "Brain Gene" = "hotpink",
      "Non-brain Gene" = scales::alpha("hotpink", 0.45)
    )
  ) +
  labs(
    title = "B",subtitle = "NeuN",
    y = "Peak set",
    x = "% of annotated genes",
    fill = "Gene annotation"
  ) +
  pub_theme +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12) ,#angle = 45, hjust =1),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 24,face = "bold",hjust = -0.25),
    plot.subtitle = element_text(size = 22,
                                 face = "bold",
                                 hjust = 0.5)
    
  )+ guides(fill = guide_legend(ncol = 2))

brainGenes.plot.SOX <- brainGenes.plot[brainGenes.plot$cellType == "SOX10",]
level_order <- c("Dataset 1 specific default", "Dataset 1 shared default",
                 "Dataset 1 specific curated","Dataset 1 shared curated",
                 "Dataset 2 specific default", "Dataset 2 shared default",
                 "Dataset 2 specific curated","Dataset 2 shared curated")
brainGenes.plot.SOX$peakset <- factor(brainGenes.plot.SOX$peakset,
                                      levels =rev(level_order ))
p14B2 <- ggplot(
  brainGenes.plot.SOX,
  aes(
    y = peakset,
    x = num,
    fill = annotation
  )
) +
  geom_col(
    position = "fill",
    width = 0.8,
    linewidth = 0,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) x * 100
  )+ 
  scale_fill_manual(
    values = c(
      "Brain Gene" = "#6B8E23",
      "Non-brain Gene" = scales::alpha("#6B8E23", 0.45)
    )
  ) +
  labs(
    title = " ",subtitle = "SOX10",
    y = "Peak set",
    x = "% of annotated genes",
    fill = "Gene annotation"
  ) +
  pub_theme +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12) ,#angle = 45, hjust =1),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 24,face = "bold",hjust = -0.25),
    plot.subtitle = element_text(size = 22,
                                 face = "bold",
                                 hjust = 0.5)
    
  )+ guides(fill = guide_legend(ncol = 2))

panel_top <- (p14A1 | p14A2) +
  plot_layout(widths = c(1.5,1.5))

middle_pannel <- p14B1 | p14B2

SupFig14 <- panel_top / middle_pannel  +
  plot_layout(heights = c(0.75,1))


ggsave(
  filename = "SupplementaryFig_14.png",
  plot = SupFig14,
  width = 20,
  height = 15
)

ggsave(
  filename = "SupplementaryFig_14.pdf",
  plot = SupFig14,
  width = 20,
  height = 15
)
                         
