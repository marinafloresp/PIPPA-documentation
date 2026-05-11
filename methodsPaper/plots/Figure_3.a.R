### Figure 3a
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures")
#Line graph of number of peaks (y-axis) vs q-value threshold (x-axis). BDR & SFARI, three positive cell types, one line per dataset (6 in total). 
library(GenomicRanges)
#Dataset 1 Q5e-2 peaks
neun_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/NEUN_5e-2.sorted.narrowPeak.filt")
sox10_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/SOX_5e-2.sorted.narrowPeak.filt")
irf8_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/IRF8_5e-2.sorted.narrowPeak.filt")
##Dataset 2 Q5e-2 peaks
neun_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/NEUN/neun_q5e-2.narrowPeak.filt")
sox10_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/SOX10/sox10_q5e-2.narrowPeak.filt")
irf8_dat2 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/IRF8/irf8_q5e-2.narrowPeak.filt")



n_vals <- c(-log10(0.05), 2:100)
#Function to count peaks passing threshold

count_peaks <- function(gr, n_vals) {
  # If qValue exists, count directly
  if("qValue" %in% colnames(mcols(gr))) {
    sapply(n_vals, function(n) sum(mcols(gr)$qValue >= n, na.rm = TRUE))
  } else {
    stop("GRanges has no qValue column!")
  }
}

# count per threshold per dataset
neun_dat1_counts     <- count_peaks(neun_dat1, n_vals)
sox10_dat1_counts   <- count_peaks(sox10_dat1, n_vals)
irf8_dat1_counts <- count_peaks(irf8_dat1, n_vals)
neun_dat2_counts     <- count_peaks(neun_dat2, n_vals)
sox10_dat2_counts   <- count_peaks(sox10_dat2, n_vals)
irf8_dat2_counts <- count_peaks(irf8_dat2, n_vals)

#create dataframe of peak theshold counts
df <- data.frame(
  n = rep(n_vals, 6),
  peaks = c(
    neun_dat1_counts,
    neun_dat2_counts,
    sox10_dat1_counts,
    sox10_dat2_counts,
    irf8_dat1_counts,
    irf8_dat2_counts
  ),
  dataset = rep(
    c(
      "NeuN Datastet.1",
      "NeuN Dataset.2",
      "SOX10 Dataset.1",
      "SOX10 Dataset.2",
      "IRF8 Dataset.1",
      "IRF8 Dataset.2"
    ),
    each = length(n_vals)
  )
)
head(df)

custom_colors <- c(
  "NeuN Datastet.1" = "deeppink4",        
  "SOX10 Dataset.1" = "darkslategrey", 
  "IRF8 Dataset.1" = "darkblue",
  "NeuN Dataset.2" = "indianred",
  "SOX10 Dataset.2" = "darkslategray4",#
  "IRF8 Dataset.2" = "cornflowerblue"
)

# -----------------------------
# 7. Plot
# -----------------------------



p <- ggplot(df, aes(x = n, y = peaks, color = dataset)) +
  geom_line(size = 0.8) +
  geom_point(size = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = seq(2, 100, 4),
    labels = function(x) paste0("1e-", x)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Q Value",
    y = "Number of Peaks",
    title = "Peak counts vs Q value threshold per dataset",
    color = "Peak Set"
  )


ggsave("Figure_5.a.pdf", height = 8, width = 10)