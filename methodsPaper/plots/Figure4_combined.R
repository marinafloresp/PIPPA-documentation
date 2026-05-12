## Figure 4
#Figure 4: More stringent q-value thresholds change the genomic annotation of peaks. BDR & SFARI, three positive cell types, one line per dataset (6 in total). 

# 4a Line graph of number of genes annotated to any peak (y-axis) vs qvalue threshold. 

# 4b Line graph of number of peaks annotated to each gene (y-axis) vs qvalue threshold. 

# 4c Line graph/stacked barchart of % peaks annotated to promoter, exon, intron, intergenic vs qvalue threshold. 


# libraries
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

#Dataset 1 Q5e-2 peaks
neun_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/NEUN_5e-2.sorted.narrowPeak.filt")
sox10_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/SOX_5e-2.sorted.narrowPeak.filt")
irf8_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/IRF8_5e-2.sorted.narrowPeak.filt")
##Dataset 2 Q5e-2 peaks
neun_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/NEUN/neun_q5e-2.narrowPeak.filt")
sox10_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/SOX10/sox10_q5e-2.narrowPeak.filt")
irf8_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/IRF8/irf8_q5e-2.narrowPeak.filt")


custom_colors <- c(
  "NeuN Datastet.1" = "deeppink4",        
  "SOX10 Dataset.1" = "darkslategrey", 
  "IRF8 Dataset.1" = "darkblue",
  "NeuN Dataset.2" = "indianred",
  "SOX10 Dataset.2" = "darkslategray4",#
  "IRF8 Dataset.2" = "cornflowerblue"
)


#Figure 4a

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# thresholds as -log10(q)
n_vals <- c(-log10(0.05), 2:100)

## annotate peaks per dataset
### neun ###
peakAnno <- annotatePeak(
  neun_dat1,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df <- as.data.frame(peakAnno)

peakAnno_neun2 <- annotatePeak(
  neun_dat2,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df_neun2 <- as.data.frame(peakAnno_neun2)

### sox10 ###
peakAnno_sox10.1 <- annotatePeak(
  sox10_dat1,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df_sox10.1 <- as.data.frame(peakAnno_sox10.1)

peakAnno_sox10.2 <- annotatePeak(
  sox10_dat2,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df_sox10.2 <- as.data.frame(peakAnno_sox10.2)


### irf8 ###
peakAnno_irf8.1 <- annotatePeak(
  irf8_dat1,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df_irf8.1 <- as.data.frame(peakAnno_irf8.1)

peakAnno_irf8.2 <- annotatePeak(
  irf8_dat2,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anno_df_irf8.2 <- as.data.frame(peakAnno_irf8.2)


  
  process_thresholds <- function(anno_df, sample_name, n_vals){
    
    results <- list()
    
    for(n in n_vals){
      
      # subset by threshold
      df_sub <- anno_df %>%
        filter(qValue >= n)
      
      if(nrow(df_sub) == 0) next
      
      # 4a unique genes
      n_genes <- df_sub %>%
        mutate(SYMBOL = as.character(SYMBOL)) %>%
        filter(!is.na(SYMBOL), SYMBOL != "") %>%
        distinct(SYMBOL) %>%
        nrow()
      
      # 4b peaks per gene
      mean_peaks_per_gene <- df_sub %>%
        mutate(SYMBOL = as.character(SYMBOL)) %>%
        filter(!is.na(SYMBOL), SYMBOL != "") %>%
        group_by(SYMBOL) %>%
        summarise(n = n(), .groups = "drop") %>%
        summarise(mean_peaks = mean(n)) %>%
        pull(mean_peaks)
      
      # 4c annotation proportions
      annotation_df <- df_sub %>%
        mutate(
          annotation_group = case_when(
            grepl("Promoter", annotation) ~ "Promoter",
            grepl("Exon", annotation) ~ "Exon",
            grepl("Intron", annotation) ~ "Intron",
            TRUE ~ "Intergenic"
          )
        ) %>%
        group_by(annotation_group) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(percent = n / sum(n) * 100)
      
      results[[as.character(n)]] <- list(
        genes = data.frame(
          sample = sample_name,
          threshold = n,
          n_genes = n_genes
        ),
        
        peaks_per_gene = data.frame(
          sample = sample_name,
          threshold = n,
          mean_peaks_per_gene = mean_peaks_per_gene
        ),
        
        annotation = annotation_df %>%
          mutate(
            sample = sample_name,
            threshold = n
          )
      )
    }
    
    return(results)
  }
## Generate results per peak set
# neun, D1
  neun1_results <- process_thresholds(
    anno_df,
    "NeuN Dataset.1",
    n_vals
  )
  
# neun D2
  neun2_results <- process_thresholds(
    anno_df_neun2,
    "NeuN Dataset.2",
    n_vals
  )
  
  # sox10 D1
  sox10.1_results <- process_thresholds(
    anno_df_sox10.1,
    "SOX10 Dataset.1",
    n_vals
  )
  
  # sox10 D2
  sox10.2_results <- process_thresholds(
    anno_df_sox10.1,
    "SOX10 Dataset.2",
    n_vals
  )
  
  # irf8 D1
  irf8.1_results <- process_thresholds(
    anno_df_irf8.1,
    "IRF8 Dataset.1",
    n_vals
  )
  
  # irf8 D2
  irf8.1_results <- process_thresholds(
    anno_df_irf8.2,
    "IRF8 Dataset.2",
    n_vals
  )
  
summary(neun1_results)

# Fig 4a
fig4a_df <- bind_rows(
  lapply(neun1_results, function(x) x$genes),
  lapply(neun2_results, function(x) x$genes),
  lapply(irf8.1_results, function(x) x$genes)
)

ggplot(fig4a_df,
       aes(x = threshold,
           y = n_genes,
           color = sample,
           group = sample)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = "q-value threshold",
    y = "Number of annotated genes"
  )


## Fig 4b
fig4b_df <- bind_rows(
  results_neun1$peaks_per_gene
  # add all others here
)

ggplot(fig4b_df,
       aes(x = threshold,
           y = mean_peaks_per_gene,
           color = sample,
           group = sample)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = "q-value threshold",
    y = "Mean peaks per gene"
  )

## Fig 4C
fig4c_df <- bind_rows(
  results_neun1$annotation
  # add all others here
)

ggplot(fig4c_df,
       aes(x = threshold,
           y = percent,
           fill = annotation_group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~sample) +
  theme_classic() +
  labs(
    x = "q-value threshold",
    y = "% Peaks",
    fill = "Annotation"
  )


## updated
fig4a_df <- bind_rows(
  lapply(neun1_results, function(x) x$genes)
)

fig4b_df <- bind_rows(
  lapply(neun1_results, function(x) x$peaks_per_gene)
)

fig4c_df <- bind_rows(
  lapply(neun1_results, function(x) x$annotation)
)