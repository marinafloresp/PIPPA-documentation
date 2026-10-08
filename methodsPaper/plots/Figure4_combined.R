#########################
##       Figure 4      ##
#########################

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
  "NeuN Dataset.1" = "deeppink4",        
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
  
  for(thresh in n_vals){

    # subset by threshold
    df_sub <- anno_df %>%
      dplyr::filter(qValue >= thresh)
    
    if(nrow(df_sub) == 0) next
    # 4A number of unique genes
    n_genes <- df_sub %>%
      dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
      dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
      dplyr::distinct(SYMBOL) %>%
      nrow()
    # 4B peaks per gene
    
    mean_peaks_per_gene <- df_sub %>%
      dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
      dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(
        peak_count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::summarise(
        mean_peaks = mean(peak_count)
      ) %>%
      dplyr::pull(mean_peaks)
    
    # 4C annotation proportions
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
    

    # store
    
    results[[as.character(thresh)]] <- list(
      
      genes = data.frame(
        Peaks = sample_name,
        Threshold = thresh,
        n_genes = n_genes
      ),
      
      peaks_per_gene = data.frame(
        Peaks = sample_name,
        Threshold = thresh,
        mean_peaks_per_gene = mean_peaks_per_gene
      ),
      
      annotation = annotation_df %>%
        dplyr::mutate(
          Peaks = sample_name,
          Threshold = thresh
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
    anno_df_sox10.2,
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
  irf8.2_results <- process_thresholds(
    anno_df_irf8.2,
    "IRF8 Dataset.2",
    n_vals
  )
  
#summary(neun1_results)

  all_results <- list(
    "NeuN.1" = neun1_results,
    "NeuN.2" = neun2_results,
    "SOX10.1" = sox10.1_results,
    "SOX10.2" = sox10.2_results,
    "IRF8.1" = irf8.1_results,
    "IRF8.2" = irf8.2_results
  )
# Fig 4a
fig4a_df <- bind_rows(
  lapply(all_results, function(res_list){
    bind_rows(lapply(res_list, `[[`, "genes"))
  })
)

head(fig4a_df)
table(fig4a_df$Peaks)
fig4a_df$DatasetType <- ifelse(
  grepl("Dataset.2", fig4a_df$Peaks),
  "Dataset.2",
  "Dataset.1"
)
colnames(fig4a_df) <- c("Peaks", "Threshold", "n_genes",  "Dataset")
### plot fig4a
library(ggplot2)
fig4a <- ggplot(fig4a_df,
                aes(x = Threshold,
                    y = n_genes,
                    color = Peaks,
                    linetype = Dataset,
                    group = Peaks)) +
  geom_line(linewidth = 0.7) +
  #geom_point(size = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c(
    "Dataset.1" = "solid",
    "Dataset.2" = "dashed"
  )) +
  theme_classic() +
  labs(
    #x = "Q value, 1e-n",
    x = expression(-log[10]~"(Q value)"),
    y = "Total annotated genes"
  )

ggsave("Figure_4a.pdf", height = 7, width = 7)

# plot
fig4a_df <- as.data.frame(fig4a_df)
head(fig4a_df)
# 4a Line graph of number of genes annotated to any peak (y-axis) vs qvalue threshold. 


########################## Fig 4b ###########################################

fig4b_df <- bind_rows(
  lapply(all_results, function(res_list){
    bind_rows(lapply(res_list, `[[`, "peaks_per_gene"))
  })
)

fig4b_df <- as.data.frame(fig4b_df)

fig4b_df$Dataset <- ifelse(
  grepl("Dataset.2", fig4b_df$sample),
  "Dataset.2",
  "Dataset.1"
)


colnames(fig4b_df) <- c("Peaks", "Threshold", "mean_peaks_per_gene",  "Dataset")

head(fig4b_df)

fig4b <- ggplot(fig4b_df,
       aes(x = Threshold,
           y = mean_peaks_per_gene,
           color = Peaks,
           group = Peaks,
           linetype = Dataset)) +
  geom_line(linewidth = 0.7) +
#  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c(
    "Dataset.1" = "solid",
    "Dataset.2" = "dashed"
  )) +
  theme_classic() +
  labs(
    #x = "Q value, 1e-n",
    x = expression(-log[10]~"(Q value)"),
    y = "Mean peaks per gene"
  )

ggsave("Figure_4b.pdf", height = 7, width = 7)

## Fig 4C
fig4c_df <- bind_rows(
  lapply(all_results, function(res_list){
    bind_rows(lapply(res_list, `[[`, "annotation"))
  })
)
fig4c_df <- as.data.frame(fig4c_df)
head(fig4c_df)

n_vals <- c(-log10(0.05), seq(10, 100, by = 10))
# rerun process_function with n_vals updated
fig4c <- ggplot(fig4c_df,
       aes(x = Threshold,
           y = percent,
           fill = annotation_group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Peaks) +
  theme_classic() +
  labs(
    x = expression(-log[10]~"(Q value)"),
    y = "% Peaks",
    fill = "Annotation"
  )

ggsave('Figure_4c.pdf', width = 14, height = 7)


### panel plot

#  combined Figure 4 panel


library(patchwork)
library(cowplot)


# Shared theme


pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )



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




# FIGURE 4A


fig4a <- ggplot(
  fig4a_df,
  aes(
    x = Threshold,
    y = n_genes,
    color = Peaks,
    linetype = Dataset,
    group = Peaks
  )
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c(
    "Dataset.1" = "solid",
    "Dataset.2" = "dashed"
  )) +
  scale_x_continuous(
    breaks = c(-log10(0.05), seq(10, 100, by = 10)),
    labels = c("1.3", seq(10, 100, by = 10))
  ) +
  labs(
    x = expression(-log[10]~"(Q value)"),
    y = "Annotated genes",
    title = "A"
  ) +
  pub_theme


# FIGURE 4B


fig4b <- ggplot(
  fig4b_df,
  aes(
    x = Threshold,
    y = mean_peaks_per_gene,
    color = Peaks,
    linetype = Dataset,
    group = Peaks
  )
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c(
    "Dataset.1" = "solid",
    "Dataset.2" = "dashed"
  )) +
  scale_x_continuous(
    breaks = c(-log10(0.05), seq(10, 100, by = 10)),
    labels = c("1.3", seq(10, 100, by = 10))
  ) +
  labs(
    x = expression(-log[10]~"(Q value)"),
    y = "Mean peaks per gene",
    title = "B"
  ) +
  pub_theme

# FIGURE 4C


fig4c <- ggplot(
  fig4c_df,
  aes(
    x = factor(Threshold),
    y = percent,
    fill = annotation_group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  facet_wrap(~Peaks, nrow = 2) +
  scale_x_discrete(
    labels = c("1.3", seq(10, 100, by = 10))
  ) +
  labs(
    x = expression(-log[10]~"(Q value)"),
    y = "% Peaks",
    fill = "Genomic annotation",
    title = "C"
  ) +
  pub_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# Extract shared legend from fig4a


legend_shared <- cowplot::get_legend(
  fig4a +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
)

# Remove legends from individual panels


fig4a_clean <- fig4a + theme(legend.position = "none")
fig4b_clean <- fig4b + theme(legend.position = "none")
fig4c_clean <- fig4c #+ theme(legend.position = "none")


# Combine top row


top_row <- fig4a_clean | fig4b_clean


# Final combined figure


combined_fig4 <-
  (top_row / fig4c_clean / legend_shared) +
  plot_layout(
    heights = c(1, 1.4, 0.15)
  )


# Save publication figure


ggsave(
  filename = "Figure4_combined.pdf",
  plot = combined_fig4,
  width = 14,
  height = 16
)

ggsave(
  filename = "Figure4_combined.png",
  plot = combined_fig4,
  width = 14,
  height = 16,
  dpi = 600
)




