####################################
##       Supplementary Fig 9      ##
####################################

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(forcats)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(ensembldb)

## brain genes list ##
brainGenes <- read.table("./brainGenes.txt")[,1]
brainGenes <- unique(brainGenes)
brainGenes <- sub("\\..*$", "", brainGenes)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb <- EnsDb.Hsapiens.v86

groupsPeaks <-"./methodsPaper/data/downReads/"

## plotting settings ##
pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 14),
    legend.key.width = unit(2.2, "cm")
  )

col_vals <- c("0.05" = "#56B4E9","0.01" = "salmon","0.00001" = "#009E73")
shape_vals_ct <- c("NeuN+" = "solid", "SOX10+" = "dotted")

## Get protein coding genes IDs ##
protein_coding_ids <- genes(edb, 
                             filter = GeneBiotypeFilter("protein_coding") )$gene_id 
brainGenes.protein <- brainGenes[brainGenes %in% protein_coding_ids]

######################
## Gather peak sets ##
######################

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


calculatePeaksGenes <- function(peakAnno, setName) { 
  annot <- peakAnno %>% as.data.frame() 
  peaks_per_gene <- annot %>% 
    dplyr::filter( !is.na(SYMBOL), SYMBOL != "", 
                   ENSEMBL %in% protein_coding_ids ) %>% 
    dplyr::count(SYMBOL, name = "peak_count") %>% 
    dplyr::mutate(peak_set = setName) 
  list( summary = tibble::tibble( name = setName, 
                                  n_genes = nrow(peaks_per_gene), 
  mean_peaks_per_gene = mean(peaks_per_gene$peak_count) ), 
  peaks_per_gene = peaks_per_gene ) 
  }

results <- purrr::imap( peakAnno.all, ~ calculatePeaksGenes( peakAnno = .x, setName = .y) )
all_summaries <- lapply(results, `[[`, "summary")

genesReadsSets.all <- bind_rows(all_summaries)
genesReadsSets.all$num_reads <- as.integer(str_extract(genesReadsSets.all$name, "(?<=_)[0-9]+(?=readsAll)"))
genesReadsSets.all$thres <- str_extract(genesReadsSets.all$name, "[0-9]+e-[0-9]+")
genesReadsSets.all$fraction <-str_extract(genesReadsSets.all$name, "^[^_]+")
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "NEUN"] <- "NeuN+"
genesReadsSets.all$fraction[genesReadsSets.all$fraction == "SOX"] <- "SOX10+"
genesReadsSets.all$thres <-fct_recode(genesReadsSets.all$thres, "0.05" = "5e-2", "0.01" = "1e-2", "0.00001" = "1e-5")

# ########################## Fig Sup 9.a ###########################################
# Line plot of total annotated genes (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold
# line type by cell type

p9a <- ggplot(as.data.frame(genesReadsSets.all), aes(x = num_reads, y = n_genes, color = thres, linetype = fraction)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSets.all$num_reads))
  ) +
  scale_linetype_manual(values = shape_vals_ct) + 
  scale_color_manual(values = col_vals)+
  labs(
    x =expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Total annotated genes",
    title = "A",
    color = "Q value threshold", linetype = "Cell type"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 9.b ###########################################
# Line plot of mean genes annotated per peak (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold
# line type by cell type

p9b <- ggplot(
  genesReadsSets.all,
  aes(x = factor(num_reads), y = mean_peaks_per_gene, color = thres,
      linetype = fraction,
      group = interaction(thres, fraction))
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_discrete(
    breaks = sort(unique(genesReadsSets.all$num_reads))
  )+
  # scale_y_log10()+
  scale_linetype_manual(values = shape_vals_ct) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Mean peaks per gene", 
    color = "Q value threshold", 
    linetype = "Cell type",
    title = "B"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 9.c ###########################################
# Line plot of % of peaks annotated as brain genes (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold
# line type by cell type

calculateBrainGenes <- function(peakAnno, setName){
  annot <- peakAnno@anno %>%
    as.data.frame()
  
  annot$ENSEMBL <- sub("\\..*$", "", annot$ENSEMBL)
  
  ids <- unique(annot$ENSEMBL)
  genes_with_peaks <- ids
  intersect <- intersect(genes_with_peaks,brainGenes.protein)
  
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

p9c <- ggplot(
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
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = sort(unique(genesReadsSetsBrain.all$num_reads))
  ) +
  scale_linetype_manual(values = shape_vals_ct) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "% Peaks annotated as brain genes",
    color = "Q value threshold", linetype = "Cell type",
    title = "C"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

# ########################## Fig Sup 9.d ###########################################
# Line plot of mean number of narrow (< 400bp) peaks annotated per gene (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold
# line type by cell type

calculatePeaksGenes <- function(peakAnno, setName){
  
  annot <- peakAnno %>%
    as.data.frame() %>%
    dplyr::mutate(
      peak_size = end - start + 1,
      short_peak = peak_size < 400,
      promoter_peak = grepl("Promoter", annotation)
    )
  peaks_per_gene <- annot %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "", 
                  ENSEMBL %in% protein_coding_ids) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      peak_count = dplyr::n(),
      peaks_lt400 = sum(short_peak),
      promoter_peaks = sum(promoter_peak),
      fraction_lt400 = peaks_lt400 / peak_count,
      fraction_promoter = promoter_peaks / peak_count,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      peak_set = setName
    
  list(
    summary = tibble::tibble(
      name = setName,
      n_genes = nrow(peaks_per_gene),
      mean_peaks_per_gene = mean(peaks_per_gene$peak_count),
      total_annotated_peaks = sum(peaks_per_gene$peak_count),
      total_peaks_lt400 = sum(peaks_per_gene$peaks_lt400),
      fraction_peaks_lt400 =
        sum(peaks_per_gene$peaks_lt400) /
        sum(peaks_per_gene$peak_count),
      mean_peaks_lt400 =
        mean(peaks_per_gene$peaks_lt400),
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

p9d <- ggplot(
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
  scale_linetype_manual(values = shape_vals_ct) + 
  scale_color_manual(values = col_vals)+
  labs(
    x =expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Mean narrow peaks (<400bp) per gene",
    color = "Q value threshold", linetype = "Cell type",
    title = "D"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )
      
# ########################## Fig Sup 9.e ###########################################
# Line plot of mean number of promoert peaks annotated per gene (y-axis)
# by read pairs (x-axis) 
# coloured by Q value threshold
# line type by cell type
      
p9e <- ggplot(
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
  scale_linetype_manual(values = shape_vals_ct) + 
  scale_color_manual(values = col_vals)+
  labs(
    x = expression("Read pairs (" ~ x ~ 10^7 ~ ")"),
    y = "Mean promoter peaks per gene",
    color = "Q value threshold", linetype = "Cell type",
    title = "E"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom"
  )

## combine plots into figure ##
legend_shared <- cowplot::get_legend(
  p9a +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
)
p9a_clean <- p9a + theme(legend.position = "none")
p9b_clean <- p9b2 + theme(legend.position = "none")

p9c_clean <- p9c + theme(legend.position = "none")
p9d_clean <- p9d + theme(legend.position = "none")
p9e_clean <- p9e + theme(legend.position = "none")
legend_shared <- cowplot::get_legend(
  p9a +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
)

# Combine top row

figSup9 <- 
  (p9a_clean | p9b_clean) / (p9c_clean | p9d_clean) / (p9e_clean | legend_shared) +
  plot_layout(
    heights = c(1, 1, 1), nrow = 3
  )

ggsave(
  filename = "SupplementaryFig_9_D1_proteinCoding.pdf",
  plot = figSup9,
  width = 14,
  height = 16
)

ggsave(
  filename = "SupplementaryFig_9_D1_proteinCoding.png",
  plot = figSup9,
  width = 14,
  height = 16
)

