#####################################
##       Supplementary Fig 4       ##
#####################################

options(scipen=5)
library(ggplot2)
library(ATACseqQC)
library(plyr)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

## settings for plotting ##

pub_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11 ),
    plot.title = element_text(face = "bold", size = 24),
    legend.position = "bottom",
    legend.key.width = unit(2.2, "cm")
  )

###############
## functions ##
###############

fdsPlots <- function(listSamples){
  counter = 1
  list_plots <-list()
  samples <- names(listSamples)
  for(i in 1:length(listSamples)) {
    df <-data.frame(listSamples[i])
    colnames(df)<- c("Var1", "Freq")
    df$sample <-samples[i]
    list_plots[[i]] <- df
    counter = counter + 1
  }
  bigdf <-ldply(list_plots, data.frame)
  bigdf[,1]<- as.integer(as.character(bigdf[,1]))
  return(bigdf)
}

samples2Plot <- read.table("list_of_samples.txt")[,1]
qcStats <- read.csv(file.path(qc1stats),header = TRUE, sep = ',', stringsAsFactors = FALSE)
qcMetrics <- qcMetrics[match(samples2Plot,qcMetrics$sampleID), ]
qcMetrics$QC_status <- rep(c("PASS", "FAIL"), each = 3)

samples <- samples2Plot
aQCFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam.bai$", recursive = TRUE, full.names = TRUE)
aQCFiles<-gsub(".bai", "", aQCFiles)
aQCFiles <- aQCFiles[match(samples,gsub(paste0(alignedDir,"/"),"",gsub(".filt.nodup.bam", "", aQCFiles)))]
aQCSampleNames<-gsub(".filt.nodup.bam", "", basename(aQCFiles))
aQCSampleNames <- na.omit(aQCSampleNames)
aQCFiles <- na.omit(aQCFiles)
nSamples <- length(samples)

# ########################## Fig Sup 4 ############################
# Fragment size distribution of pass vs failed samples 

if(nSamples > 0){

	## create summary of fragment size using filtered aligned files
  print("Samples for fragment size distribution: ")
  print(samples)
	fragSizeHist<-fragSizeDist(aQCFiles, samples)
  
  save(fragSizeHist, file = "FragSize_periodicity_vals.rdata")
  load("FragSize_periodicity_vals.rdata")
  BatchSamples <- fdsPlots(fragSizeHist)
  
  BatchSamples$sample <- factor(
  BatchSamples$sample,
  levels = samples2Plot
  )
  BatchSamples.split <- split(BatchSamples,BatchSamples$sample)
  BatchSamples.split <- split(
    BatchSamples,
    BatchSamples$sample
  )

  p<- NULL
  
  for(i in 1:length(BatchSamples.split)){
  
  fdValues <- BatchSamples.split[[i]]
  axis1 <- seq(from = 0, to = max(fdValues[,1]), by = 100)
  
  qc_text <- paste0(
    "NFR: ", round(qcMetrics$NFR[i],3),
    "\nMNC: ", round(qcMetrics$MNC[i],3)
  )
  
  p[[i]] <- ggplot(fdValues, aes(x = Var1, y = Freq)) +
    geom_line(linewidth =0.5) +
    pub_theme +
    labs(
      x = "Fragment size (bp)",
      y = "Frequency"
    ) +
    scale_x_continuous(
      breaks = axis1,
      limits = c(0, 600)
    ) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = qc_text,
      hjust = 1.1,
      vjust = 1.1,
      size = 5,
      fill = ifelse(qcMetrics$QC_status[i] == "PASS",
                    "#009E73",
                    "#D55E00")
    )
  }
  
  pdf(file.path(paste0("SupplementaryFig_4_D1.pdf")), width = 20, height = 10)
  grid.arrange(grobs = p, ncol=3)
  dev.off()
  
  png(file.path(paste0("SupplementaryFig_4_D1.png")), width = 2000, height = 1000)
   grid.arrange(grobs = p, ncol=3)
  dev.off()
}
