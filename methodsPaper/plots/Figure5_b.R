

###############################################################################################
############################### segmented analysis thresholds #################################
###############################################################################################
## plot peak numbers:
library(GenomicRanges)
library(ggplot2)

## running cell types separately for the first step using peak set names in order to fix the different name syntax issues
# NeuN
## Dataset 1
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/BDR_threshPeaks/segmented")
peak_files <- list.files(
  path = ".",
  pattern = "^NEUN.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)

names(peak_list) <- paste0("Dat1_", names(peak_list))

names(peak_list) <- gsub(
  ".sorted",
  "",
 names(peak_list)
)

names(peak_list) <- gsub(
  ".pva",
  "",
  names(peak_list)
)

peak_list_n1 <- peak_list

# Inspect names
names(peak_list)
## Dataset 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/SFARI_threshPeaks/segmented")
# Find files
peak_files <- list.files(
  path = ".",
  pattern = "^neun.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)
names(peak_list) <- paste0("Dat2_", names(peak_list))
# Inspect
names(peak_list)

peak_list_n2 <- peak_list

## combine lists
peak_list<- c(peak_list_n1, peak_list_n2)
names(peak_list)
names(peak_list) <- gsub("NEUN", "NeuN", names(peak_list))
names(peak_list) <- gsub("neun", "NeuN", names(peak_list))
peak_list_neun <- peak_list

peak_sets <- peak_list_neun

df <- data.frame(
  dataset = names(peak_sets),
  n_peaks = sapply(peak_sets, length),
  min_log10_q = sapply(peak_sets, function(gr) {
    min(mcols(gr)$qValue, na.rm = TRUE)
  })
)

# Convert back to raw q-value threshold
df$q_thresh <- 10^(-df$min_log10_q)


# Clean labels for plotting

df$label <- paste0(
  "q<", format(df$q_thresh, scientific = TRUE, digits = 2),
  "\n(-log10=", round(df$min_log10_q, 1), ")"
)

df$dataset <- gsub("NEUN","NeuN", df$dataset)
df$dataset <- gsub("neun","NeuN", df$dataset)
df$dataset <- gsub("1e","Q1e", df$dataset)
df$dataset <- gsub("q1e","Q1e", df$dataset)
df$dataset <- gsub("qQ","Q", df$dataset)
# order by peak count 
df$dataset <- factor(df$dataset, levels = df$dataset[order(df$n_peaks)])
df$dataset_group <- sub(
  "_.*",
  "",
  df$dataset
)

head(df[, c("dataset", "dataset_group")])
head(df)

# Plot
# Clean x-axis labels

df$dataset_label <- gsub(
  "_",
  " ",
  df$dataset
)

p1 <- ggplot(
  df,
  aes(
    x = dataset_label,
    y = n_peaks,
    fill = dataset_group,
  )
) +
  
  geom_col(
    linewidth = 0.7,
    width = 0.8,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Dat1" = "hotpink",
      "Dat2" = scales::alpha("hotpink", 0.45)
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Peak Set",
    y = "Number of Peaks",
    title = "NeuN"
  ) +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    
    legend.position = "none"
  )

pdf("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/peakNumber_curves/segmented/peak_numbers_neun.pdf")
p1
dev.off()
df$cell <- "NeuN"
df_num_neun <- df
head(df_num_neun)
#########.     SOX10. ##################
## Dataset 1
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/BDR_threshPeaks/segmented")
peak_files <- list.files(
  path = ".",
  pattern = "^SOX.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)

names(peak_list) <- paste0("Dat1_", names(peak_list))

names(peak_list) <- gsub(
  ".sorted",
  "",
  names(peak_list)
)

names(peak_list) <- gsub(
  ".pva",
  "",
  names(peak_list)
)

peak_list_n1 <- peak_list

# Inspect names
names(peak_list)
## Dataset 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/SFARI_threshPeaks/segmented")
# Find files
peak_files <- list.files(
  path = ".",
  pattern = "^sox10.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)
names(peak_list) <- paste0("Dat2_", names(peak_list))
# Inspect
names(peak_list)

peak_list_n2 <- peak_list

## combine lists
peak_list <- c(peak_list_n1, peak_list_n2)
names(peak_list)
names(peak_list) <- gsub("SOX", "SOX10", names(peak_list))
names(peak_list) <- gsub("sox10", "SOX10", names(peak_list))
peak_list_sox10 <- peak_list

### plot peak numbers per cell fraction across datasets
## plot peak numbers:

# ==============================
# Extract peak counts + thresholds
# ==============================
peak_sets <- peak_list_sox10
df <- data.frame(
  dataset = names(peak_sets),
  n_peaks = sapply(peak_sets, length),
  min_log10_q = sapply(peak_sets, function(gr) {
    min(mcols(gr)$qValue, na.rm = TRUE)
  })
)


head(df)

# Convert back to raw q-value threshold
df$q_thresh <- 10^(-df$min_log10_q)

# ==============================
# Clean labels for plotting
# ==============================
df$label <- paste0(
  "q<", format(df$q_thresh, scientific = TRUE, digits = 2),
  "\n(-log10=", round(df$min_log10_q, 1), ")"
)

# Optional: order by peak count (nicer plot)
df$dataset <- factor(df$dataset, levels = df$dataset[order(df$n_peaks)])
head(df)
tail(df)

#df$dataset <- gsub("SOX","SOX10", df$dataset)
df$dataset <- gsub("sox10","SOX10", df$dataset)
df$dataset <- gsub("1e","Q1e", df$dataset)
df$dataset <- gsub("q1e","Q1e", df$dataset)
df$dataset <- gsub("qQ","Q", df$dataset)

df$dataset_group <- sub(
  "_.*",
  "",
  df$dataset
)

head(df)

# Plot

# Clean x-axis labels

df$dataset_label <- gsub(
  "_",
  " ",
  df$dataset
)


## updated plot with colour contrasts
p2 <- ggplot(
  df,
  aes(
    x = dataset_label,
    y = n_peaks,
    fill = dataset_group,
  )
) +
  
  geom_col(
    linewidth = 0.7,
    width = 0.8,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Dat1" = "#6B8E23",
      "Dat2" = scales::alpha("#6B8E23", 0.45)
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Peak Set",
    y = "Number of Peaks",
    title = "SOX10"
  ) +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    
    legend.position = "none"
  )

pdf("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/peakNumber_curves/segmented/peak_numbers_sox10.pdf")
p2
dev.off()







############################    IRF8    ###################################################
## Dataset 1
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/BDR_threshPeaks/segmented")
peak_files <- list.files(
  path = ".",
  pattern = "IRF8.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)

names(peak_list) <- paste0("Dat1_", names(peak_list))

names(peak_list) <- gsub(
  ".sorted",
  "",
  names(peak_list)
)

names(peak_list) <- gsub(
  ".pva",
  "",
  names(peak_list)
)

peak_list_n1 <- peak_list

# Inspect names
names(peak_list)
## Dataset 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/SFARI_threshPeaks/segmented")
# Find files
peak_files <- list.files(
  path = ".",
  pattern = "^irf8.*\\.narrowPeak\\.filt$",
  full.names = TRUE
)

# Read all files using custom function
peak_list <- lapply(
  peak_files,
  read_narrowpeak_to_granges_signals
)

# Name list using filenames without extension
names(peak_list) <- gsub(
  "\\.narrowPeak\\.filt$",
  "",
  basename(peak_files)
)
names(peak_list) <- paste0("Dat2_", names(peak_list))
# Inspect
names(peak_list)

peak_list_n2 <- peak_list

## combine lists
peak_list <- c(peak_list_n1, peak_list_n2)
names(peak_list) <- gsub("irf8", "IRF8", names(peak_list))
names(peak_list)
summary(peak_list)
peak_list_irf8 <- peak_list

### plot peak numbers per cell fraction across datasets
## plot peak numbers:
library(GenomicRanges)
library(ggplot2)

# Extract peak counts + thresholds

peak_sets <- peak_list_irf8
df <- data.frame(
  dataset = names(peak_sets),
  n_peaks = sapply(peak_sets, length),
  min_log10_q = sapply(peak_sets, function(gr) {
    min(mcols(gr)$qValue, na.rm = TRUE)
  })
)


head(df)

# Convert back to raw q-value threshold
df$q_thresh <- 10^(-df$min_log10_q)

# Clean labels for plotting

df$label <- paste0(
  "q<", format(df$q_thresh, scientific = TRUE, digits = 2),
  "\n(-log10=", round(df$min_log10_q, 1), ")"
)

# Optional: order by peak count (nicer plot)
df$dataset <- factor(df$dataset, levels = df$dataset[order(df$n_peaks)])
head(df)

df$dataset <- gsub("irf8","IRF8", df$dataset)


##
# Clean labels for plotting
# ==============================
df$label <- paste0(
  "q<", format(df$q_thresh, scientific = TRUE, digits = 2),
  "\n(-log10=", round(df$min_log10_q, 1), ")"
)

# Optional: order by peak count (nicer plot)
df$dataset <- factor(df$dataset, levels = df$dataset[order(df$n_peaks)])
head(df)
tail(df)
df$dataset <- gsub("1e","Q1e", df$dataset)
df$dataset <- gsub("q1e","Q1e", df$dataset)
df$dataset <- gsub("qQ","Q", df$dataset)

df$dataset_group <- sub(
  "_.*",
  "",
  df$dataset
)

head(df)

# Plot

# Clean x-axis labels

df$dataset_label <- gsub(
  "_",
  " ",
  df$dataset
)

p3 <- ggplot(
  df,
  aes(
    x = dataset_label,
    y = n_peaks,
    fill = dataset_group,
  )
) +
  
  geom_col(
    linewidth = 0.7,
    width = 0.8,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Dat1" = "orange",
      "Dat2" = scales::alpha("orange", 0.45)
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Peak Set",
    y = "Number of Peaks",
    title = "IRF8"
  ) +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    
    legend.position = "none"
  )

pdf("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/peakNumber_curves/segmented/peak_numbers_irf8.pdf")
p3
dev.off()


## Combine p1, p2, p3 into a single row, make individual plot titles bold


p1 <- p1 +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p2 <- p2 +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

p3 <- p3 +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )


## Combine plots into single row for panel A


library(patchwork)

panel_A <- (p1 | p2 | p3) +
  
  plot_annotation(
    tag_levels = NULL
  ) &
  
  theme(
    plot.tag = element_text(
      size = 22,
      face = "bold"
    )
  )



############################################################
## Save
############################################################

pdf(
  "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/peakNumber_curves/segmented/panel_B_peak_numbers.pdf",
  width = 18,
  height = 6
)

panel_B

dev.off()


## add panel A and panel B
library(patchwork)
############################################################
## REMOVE OLD PATCHWORK ANNOTATIONS
############################################################

# Recreate panel_B from scratch WITHOUT plot_annotation()

panel_B <- (p1 | p2 | p3)

saveRDS(panel_B, "Fig.5.panel_B.rds")

############################################################
## PANEL A
############################################################

panel_A <-
  (
    plot_list[["NeuN_Dataset1"]] +
      plot_list[["NeuN_Dataset2"]]
  ) /
  (
    plot_list[["SOX10_Dataset1"]] +
      plot_list[["SOX10_Dataset2"]]
  ) /
  (
    plot_list[["IRF8_Dataset1"]] +
      plot_list[["IRF8_Dataset2"]]
  )

saveRDS(panel_A, "Fig.5.panel_A.rds")

############################################################
## COMBINE
############################################################

combined_plot <-
  
  wrap_elements(full = panel_A) /
  wrap_elements(full = panel_B) +
  
  plot_layout(
    heights = c(3, 1)
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
  "Figure5_combined_AB.pdf"
)

ggsave(
  filename = combined_outfile,
  plot = combined_plot,
  width = 18,
  height = 18
)





