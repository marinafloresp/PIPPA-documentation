### Figure 5a
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures")
#Line graph of number of peaks (y-axis) vs q-value threshold (x-axis). BDR & SFARI, three positive cell types, one line per dataset (6 in total). 
library(GenomicRanges)
#Dataset 1 Q5e-2 peaks
neun_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/NEUN_5e-2.sorted.narrowPeak.filt")
sox10_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/SOX_5e-2.sorted.narrowPeak.filt")
irf8_dat1 <- read_narrowpeak_to_granges_signals("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/thresholding/q0.05/BDR/200426/IRF8_5e-2.sorted.narrowPeak.filt")
##Dataset 2 Q5e-2 peaks
neun_dat2 <- 