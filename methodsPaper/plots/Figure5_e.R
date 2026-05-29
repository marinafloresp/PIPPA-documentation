## 5e annotation proportions
############################################################
## FIGURE 5E
## ChIPseeker annotation proportions
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(ChIPseeker)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno_irf8list <- lapply(peak_list_irf8 ,annotatePeak, annoDb = "org.Hs.eg.db",TxDb = txdb)
peakAnno_neunlist <-  lapply(peak_list_neun,annotatePeak, annoDb = "org.Hs.eg.db",TxDb = txdb)
peakAnno_sox10list <-  lapply(peak_list_sox10 ,annotatePeak, annoDb = "org.Hs.eg.db",TxDb = txdb)

############################################################
## FUNCTION: SUMMARISE CHIPSEEKER ANNOTATIONS
############################################################

summarise_annotations <- function(
    anno_list,
    celltype
) {
  
  bind_rows(
    
    lapply(names(anno_list), function(nm) {
      
      ######################################################
      ## EXTRACT DATAFRAME
      ######################################################
      
      df <- as.data.frame(anno_list[[nm]])
      
      ######################################################
      ## COLLAPSE ANNOTATIONS
      ######################################################
      
      df2 <- df %>%
        
        mutate(
          
          AnnotationGroup = case_when(
            
            grepl("Promoter", annotation) ~ "Promoter",
            
            grepl("Exon", annotation) ~ "Exon",
            
            grepl("Intron", annotation) ~ "Intron",
            
            TRUE ~ "Intergenic"
          )
        ) %>%
        
        group_by(AnnotationGroup) %>%
        
        summarise(
          N = n(),
          .groups = "drop"
        ) %>%
        
        mutate(
          
          Percent = N / sum(N) * 100,
          
          PeakSet = nm,
          
          CellType = celltype,
          
          Dataset = ifelse(
            grepl("^Dat1", nm, ignore.case = TRUE),
            "Dat1",
            "Dat2"
          )
        )
      
      df2
    })
  )
}

############################################################
## BUILD DATAFRAMES
############################################################

anno_neun <- summarise_annotations(
  peakAnno_neunlist,
  "NeuN"
)

anno_sox10 <- summarise_annotations(
  peakAnno_sox10list,
  "SOX10"
)

anno_irf8 <- summarise_annotations(
  peakAnno_irf8list,
  "IRF8"
)

############################################################
## CLEAN LABELS
############################################################

clean_labels <- function(x) {
  
  x <- gsub("_", " ", x)
  x <- gsub("Dat1", "Dataset 1", x)
  x <- gsub("Dat2", "Dataset 2", x)
  x <- gsub("q", "", x)
  x <- gsub("q", "", x)
  x <- gsub("1e", "Q1e", x)
  
}

anno_neun$PeakLabel  <- clean_labels(anno_neun$PeakSet)
anno_sox10$PeakLabel <- clean_labels(anno_sox10$PeakSet)
anno_irf8$PeakLabel  <- clean_labels(anno_irf8$PeakSet)

############################################################
## PLOT FUNCTION
############################################################

plot_annotation_panel <- function(
    df,
    title_text
) {
  
  ##########################################################
  ## ORDERING
  ##########################################################
  
  df$Dataset <- factor(
    df$Dataset,
    levels = c("Dat1", "Dat2")
  )
  
  df <- df %>%
   arrange(Dataset, PeakLabel)
  
  df$PeakLabel <- factor(
    df$PeakLabel,
   levels = rev(unique(df$PeakLabel))
  )
  ##########################################################
  ## PLOT
  ##########################################################
  
  ggplot(
    df,
    aes(
      x = Percent,
      y = PeakLabel,
      fill = AnnotationGroup
    )
  ) +
    
    ########################################################
  ## STACKED BARS
  ########################################################
  scale_y_discrete(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  geom_col(
    width = 0.9,
    color = NA
  ) +
    
    ########################################################
  ## COLORS
  ########################################################
  
  scale_fill_manual(
    values = c(
      "Exon" = "#F8766D",
      "Intergenic" = "#7CAE00",
      "Intron" = "#00BFC4",
      "Promoter" = "#C77CFF"
    )
  ) +
    
    ########################################################
  ## LABELS
  ########################################################
  
  labs(
    x = "% Peaks",
    y = NULL,
    title = title_text,
    fill = "Annotation"
  ) +
    
    
    ########################################################
  ## THEME
  ########################################################
  
  theme_bw(base_size = 12) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 10
      ),
      
      legend.position = "right"
    )
}

############################################################
## GENERATE PLOTS
############################################################

p_anno_neun <- plot_annotation_panel(
  anno_neun,
  "NeuN"
)

p_anno_sox10 <- plot_annotation_panel(
  anno_sox10,
  "SOX10"
)

p_anno_irf8 <- plot_annotation_panel(
  anno_irf8,
  "IRF8"
)

############################################################
## COMBINE PANEL E
############################################################

panel_E <-
  
  p_anno_neun +
  p_anno_sox10 +
  p_anno_irf8 +
  
  plot_layout(nrow = 1, guides = "collect")

############################################################
## SAVE
############################################################

ggsave(
  filename = "Figure5E_chipseeker_annotation_proportions.pdf",
  plot = panel_E,
  width = 16,
  height = 5
)

panel_E

## combine a b c d e
combined_plot.abcde <-
  
  wrap_elements(full = panel_A) /
  wrap_elements(full = panel_B) /
  wrap_elements(full=panel_C) /
  wrap_elements(full= panel_D) /
  wrap_elements(full= panel_E)+
  
  plot_layout(
    heights = c(3, 1, 1, 1, 1.3)
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag = element_text(
      size = 22,
      face = "bold"
    )
  )


############################################################
## SAVE
############################################################
outdir <-  "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/Figure_5/noElbow"
combined_outfile <- file.path(
  outdir,
  "Figure5_combined_ABCDE.pdf"
)
ggsave(
  filename = combined_outfile,
  plot = combined_plot.abcde,
  width = 18,
  height = 27
)

combined_outfile <- file.path(
  outdir,
  "Figure5_combined_ABCDE.png"
)
ggsave(
  filename = combined_outfile,
  plot = combined_plot.abcde,
  width = 18,
  height = 27
)
