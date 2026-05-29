############################################################
### Figure 5: Threshold choice analysis
### segmented regression + 95% CI
############################################################

library(kneedle)
library(segmented)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

options(scipen = 5)



peak_files <- list(
  
  "NeuN_Dataset1" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/NEUN_5e-2.sorted.narrowPeak.filt",
  
  "SOX10_Dataset1" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/SOX_5e-2.sorted.narrowPeak.filt",
  
  "IRF8_Dataset1" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/IRF8_5e-2.sorted.narrowPeak.filt",
  
  "NeuN_Dataset2" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/NEUN/neun_q5e-2.narrowPeak.filt",
  
  "SOX10_Dataset2" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/SOX10/sox10_q5e-2.narrowPeak.filt",
  
  "IRF8_Dataset2" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/IRF8/irf8_q5e-2.narrowPeak.filt"
)

peak_files

outdir <- "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/Figure_5/noElbow"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cell_colours <- c(
  "NeuN"  = "hotpink",
  "SOX10" = "#6B8E23",
  "IRF8"  = "orange"
)



## thresholds:


exponents <- -seq(2, 100, by = 2)

thresholds <- 10^exponents

############################################################
## STORE RESULTS
############################################################

results_list <- list()

summary_df <- data.frame()
plot_list <- list()



format_q_label <- function(q) {
  
  exponent <- round(log10(q))
  
  mantissa <- q / (10^exponent)
  
  nearest <- c(1, 5)[
    which.min(abs(c(1, 5) - mantissa))
  ]
  
  paste0(nearest, "e", exponent)
}



for(dataset_name in names(peak_files)) {
  
  cat("\n====================================\n")
  cat("Processing:", dataset_name, "\n")
  cat("====================================\n")
  

  ## read in peaks 

  
  peaks.gr <- makeGRangesFromDataFrame(
    read.table(peak_files[[dataset_name]]),
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    ignore.strand = TRUE
  )
  
  peaks.df <- as.data.frame(peaks.gr)
  
  colnames(peaks.df) <- c(
    "chr",
    "start",
    "end",
    "width",
    "strand",
    "peakID",
    "score",
    "strand2",
    "signal",
    "pvalue",
    "qvalue",
    "score2"
  )

  
  peaks.df$qvals <- 10^(-peaks.df$qvalue)
  

  
  counts_list <- sapply(
    thresholds,
    function(t) {
      sum(peaks.df$qvals < t, na.rm = TRUE)
    }
  )

  ## BUILD DATAFRAME

  
  df <- data.frame(
    threshold = thresholds,
    logQ = -log10(thresholds),
    count = counts_list
  )
  
  ## SEGMENTED ANALYSIS

  
  model <- lm(count ~ logQ, data = df)
  
  seg_model <- segmented(
    model,
    seg.Z = ~logQ
  )
  

  
  seg_breaks <- abs(seg_model$psi[, "Est."])
  

  
  ci_segmented <- confint(seg_model)
  
  ci_df <- as.data.frame(ci_segmented)
  
  ci_lower <- ci_df[, grep(
    "low",
    colnames(ci_df),
    ignore.case = TRUE
  )]
  
  ci_upper <- ci_df[, grep(
    "up",
    colnames(ci_df),
    ignore.case = TRUE
  )]
  

  
 
  
  seg_labels <- sapply(
    10^(-seg_breaks),
    format_q_label
  )
  
  ci_lower_labels <- sapply(
    10^(-ci_lower),
    format_q_label
  )
  
  ci_upper_labels <- sapply(
    10^(-ci_upper),
    format_q_label
  )
  

  
  results_list[[dataset_name]] <- list(
    segmented_breaks = seg_breaks,
    segmented_CI = ci_segmented
  )
  

  
  for(i in seq_len(length(seg_breaks))) {
    
    tmp <- data.frame(
      
      Dataset = dataset_name,
      
      Breakpoint_Number = i,
      
      Breakpoint_log10Q = seg_breaks[i],
      
      Breakpoint_Qvalue = signif(
        10^(-seg_breaks[i]),
        3
      ),
      
      Breakpoint_Label = seg_labels[i],
      
      CI95_Lower_log10Q = ci_lower[i],
      
      CI95_Upper_log10Q = ci_upper[i],
      
      CI95_Lower_Qvalue = signif(
        10^(-ci_upper[i]),
        3
      ),
      
      CI95_Upper_Qvalue = signif(
        10^(-ci_lower[i]),
        3
      ),
      
      CI95_Lower_Label = ci_upper_labels[i],
      
      CI95_Upper_Label = ci_lower_labels[i]
    )
    
    summary_df <- rbind(summary_df, tmp)
    
  }


  

  
  ymax <- max(df$count)
  
  seg_y <- seq(
    ymax * 0.95,
    ymax * 0.80,
    length.out = length(seg_breaks)
  )
  
  ci_low_y <- seq(
    ymax * 0.75,
    ymax * 0.65,
    length.out = length(ci_lower)
  )
  
  ci_up_y <- seq(
    ymax * 0.60,
    ymax * 0.50,
    length.out = length(ci_upper)
  )
  

  
  clean_title <- dataset_name
  
  clean_title <- gsub("_", " ", clean_title)
  
  clean_title <- gsub("Dataset1", "Dataset 1", clean_title)
  
  clean_title <- gsub("Dataset2", "Dataset 2", clean_title)
  
  # Determine colour from dataset name
  plot_colour <- case_when(
    grepl("NeuN", dataset_name, ignore.case = TRUE) ~ cell_colours["NeuN"],
    grepl("SOX10", dataset_name, ignore.case = TRUE) ~ cell_colours["SOX10"],
    grepl("IRF8", dataset_name, ignore.case = TRUE) ~ cell_colours["IRF8"]
  )
  
  
  p <- ggplot(
    df,
    aes(x = logQ, y = count)
  ) +
    
    geom_point(size = 1.3,
  color = plot_colour) +
    
    geom_line(linewidth = 1.2, color = plot_colour) +


    ########################################################
  ## Segmented breakpoints
  ########################################################
  
  
  geom_vline(
    xintercept = seg_breaks,
    color = "black",
    linewidth = 1
  ) +
    
    annotate(
      "text",
      x = seg_breaks + 0.15,
      y = ymax * 1.02,
      label = seg_labels,
      color = "black",
      vjust = 0,
      hjust = 0,
      size = 4.5,
      fontface = "bold"
    ) +
    
    ########################################################
  ## 95% CI LOWER
  ########################################################
  
  geom_vline(
    xintercept = ci_lower,
    color = "darkblue",
    linewidth = 1,
    linetype = "dashed",
    alpha = 0.7
  ) +
    
    annotate(
      "text",
      x = ci_lower - 1,
      y = ci_low_y,
      label = ci_lower_labels,
      color = "black",
      hjust = 1,
      vjust = -0.3,
      size = 4
    ) +
    
    ########################################################
  ## 95% CI UPPER
  ########################################################
  
  geom_vline(
    xintercept = ci_upper,
    color = "darkblue",
    linetype = "dashed",
    linewidth = 1,
    alpha = 0.7
  ) +
    
    annotate(
      "text",
      x = ci_upper + 1,
      y = ci_up_y,
      label = ci_upper_labels,
      color = "black",
      hjust = 0,
      vjust = -0.3,
      size = 4
    ) +
  

    
    ########################################################
  ## X AXIS
  ########################################################
  
  scale_x_continuous(
    
    name = expression(
      paste("Q value threshold (1e"^{-n}, ")")
    ),
    
    breaks = seq(2, 100, by = 2),
    
    labels = seq(2, 100, by = 2)
  ) +
    
    ########################################################
  ## LABELS
  ########################################################
  
  labs(
    y = "Number of peaks",
    title = clean_title
  ) +
    
    ########################################################
  ## THEME
  ########################################################
  
  theme_classic(base_size = 12) +
    
    theme(
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
      ),
      
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
      
    )
  
  plot_list[[dataset_name]] <- p
  
  
  ##########################################################
  ## SAVE FIGURE
  ##########################################################
  library(patchwork)
  outfile <- file.path(
    outdir,
    paste0(dataset_name, "_Figure5.pdf")
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 10,
    height = 6
  )
  
  cat("Saved figure:\n")
  cat(outfile, "\n")
}

### USE THIS ONE
# 3 by 2
combined_plot <-
  
  (
    plot_list[["NeuN_Dataset1"]] /
      plot_list[["NeuN_Dataset2"]]
  ) |
  
  (
    plot_list[["SOX10_Dataset1"]] /
      plot_list[["SOX10_Dataset2"]]
  ) |
  
  (
    plot_list[["IRF8_Dataset1"]] /
      plot_list[["IRF8_Dataset2"]]
  )

combined_outfile <- file.path(
  outdir,
  "Figure5_combined.pdf"
)

ggsave(
  filename = combined_outfile,
  plot = combined_plot,
  width = 26,
  height = 10
)

## 2 by 3


ggsave(
  filename = "combined_plot_2by3_noElbow.png",
  plot = combined_plot,
  width = 15,
  height = 14
)

saveRDS(comb_degree())
############################################################
## SAVE SUMMARY CSV
############################################################

csv_outfile <- file.path(
  outdir,
  "Figure5_threshold_summary_with95CI.csv"
)

write.csv(
  summary_df,
  file = csv_outfile,
  row.names = FALSE
)

cat("\n====================================\n")
cat("Saved summary CSV:\n")
cat(csv_outfile, "\n")
cat("====================================\n")
