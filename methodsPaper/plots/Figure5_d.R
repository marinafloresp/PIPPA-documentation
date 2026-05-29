# Figure 5 d


## FIGURE 5D
## GAP DISTANCE DISTRIBUTIONS BETWEEN PEAKS


library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)


## FUNCTION: EXTRACT INTER-PEAK GAPS


extract_peak_gaps <- function(gr_list, celltype) {
  
  bind_rows(
    
    lapply(names(gr_list), function(nm) {
      
      gr <- gr_list[[nm]]
      
      ######################################################
      ## SORT PEAKS
      ######################################################
      
      gr <- sort(gr)
      
      ######################################################
      ## CALCULATE GAPS
      ######################################################
      
      gap_gr <- gaps(gr)
      
      ######################################################
      ## KEEP STANDARD GAPS
      ######################################################
      
      gap_widths <- width(gap_gr)
      
      ######################################################
      ## REMOVE NEGATIVE / HUGE / INVALID GAPS
      ######################################################
      
      gap_widths <- gap_widths[
        is.finite(gap_widths) &
          gap_widths > 0
      ]
      
      
      ######################################################
      ## BUILD DATAFRAME
      ######################################################
      
      data.frame(
        
        PeakSet = nm,
        
        GapWidth = gap_widths,
        
        CellType = celltype,
        
        Dataset = ifelse(
          grepl("^Dat1", nm, ignore.case = TRUE),
          "Dat1",
          "Dat2"
        )
      )
    })
  )
}


## BUILD DATAFRAMES


gaps_neun <- extract_peak_gaps(
  peak_list_neun,
  "NeuN"
)

gaps_sox10 <- extract_peak_gaps(
  peak_list_sox10,
  "SOX10"
)

gaps_irf8 <- extract_peak_gaps(
  peak_list_irf8,
  "IRF8"
)


## CLEAN LABELS

clean_labels <- function(x) {
  
  x <- gsub("_", " ", x)
  x <- gsub("Dat1", "Dataset 1", x)
  x <- gsub("Dat2", "Dataset 2", x)
  x <- gsub("q", "", x)
  x <- gsub("q", "", x)
  x <- gsub("1e", "Q1e", x)
  
}

gaps_neun$PeakLabel  <- clean_labels(gaps_neun$PeakSet)
gaps_sox10$PeakLabel <- clean_labels(gaps_sox10$PeakSet)
gaps_irf8$PeakLabel  <- clean_labels(gaps_irf8$PeakSet)


## PLOT FUNCTION


plot_gap_panel <- function(
    df,
    fill_colour,
    title_text
) {
  
  ##
  ## FORCE ORDERING
  ##
  
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
  
  ##
  ## PLOT
  ##
  
  ggplot(
    df,
    aes(
      x = GapWidth,
      y = PeakLabel,
      fill = Dataset
    )
  ) +
    
    
  ## BOX PLOTS
  
  
  geom_boxplot(
    width = 0.7,
    outlier.shape = NA,
    linewidth = 1,
    color = "black"
  ) +
    
    
  ## MEDIAN POINT
  
  
  stat_summary(
    fun = median,
    geom = "point",
    size = 3,
    color = "black"
  ) +
    
    
  ## DATASET COLORS
  
  
  scale_fill_manual(
    values = c(
      "Dat1" = fill_colour,
      "Dat2" = alpha(fill_colour, 0.45)
    )
  ) +
    
    
  ## ZOOM AXIS
  
  
  coord_cartesian(
    xlim = c(
      0,
      quantile(df$GapWidth, 0.90)
    )
  ) +
    
    
  ## LABELS
  
  
  labs(
    x = "Gap between peaks (bp)",
    y = NULL,
    title = title_text
  ) +
    
    
  ## THEME
  
  
  theme_bw(base_size = 12) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      axis.text.y = element_text(
        size = 10
      ),
      
      legend.position = "none"
    )
}


## GENERATE PLOTS


p_gap_neun <- plot_gap_panel(
  gaps_neun,
  fill_colour = "hotpink",
  title_text = "NeuN"
)

p_gap_sox10 <- plot_gap_panel(
  gaps_sox10,
  fill_colour = "#6B8E23",
  title_text = "SOX10"
)

p_gap_irf8 <- plot_gap_panel(
  gaps_irf8,
  fill_colour = "orange",
  title_text = "IRF8"
)


## COMBINE PANELS


panel_D <-
  
  p_gap_neun +
  p_gap_sox10 +
  p_gap_irf8 +
  
  plot_layout(nrow = 1)


## SAVE


ggsave(
  filename = "Figure5_peak_gap_distributions.pdf",
  plot = panel_D,
  width = 15,
  height = 7
)

## combine  fig 5 A,B,C,D
combined_plot.abcd <-
  
  wrap_elements(full = panel_A) /
  wrap_elements(full = panel_B) /
  wrap_elements(full=panel_C) /
  wrap_elements(full= panel_D) +
  
  plot_layout(
    heights = c(3, 1, 1, 1)
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
  "Figure5_combined_ABCD.pdf"
)
ggsave(
  filename = combined_outfile,
  plot = combined_plot.abcd,
  width = 18,
  height = 24
)
