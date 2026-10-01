############################
##       Figure 1 b-g     ##
############################

library(edgeR)
library(ggplot2)
library(DESeq2)
library(cluster)

## universal theme for plotting ##
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

## Load metadata ##
metadata <- read.table(paste0(metaDir,"/methodsPaper/sampleSheet.csv"), header = TRUE, sep = ',', stringsAsFactors = FALSE)
dim(metadata)
metadata$sequencingBatch <- as.factor(metadata$sequencingBatch)
metadata$fraction[metadata$fraction == "SOX"] <- "Sox10+"
metadata$fraction[metadata$fraction == "NEUN"] <- "NeuN+"
metadata$fraction[metadata$fraction == "IRF8"] <- "IRF8+"
metadata$braak_stage <- metadata$clinical
metadata$clinical <- ifelse(metadata$braak_stage < 5, "control","case")
metadata$braak_stage <- as.factor(metadata$braak_stage)

## Load counts in peaks ##
load(paste0(metaDir,"/peakCounts_all_5e-2_passS12.narrowPeak.rdata"))
counts.raw <- fc_ctPeaks$counts
cnt_table <- fc_ctPeaks$annotation
colnames(counts.raw) <- sub("\\.filt\\.nodup\\.bam$", "", colnames(counts.raw))
colnames(counts.raw) <- sub("\\.sorted\\.bam$", "", colnames(counts.raw))
dim(counts.raw)

## filter metadata ##
pheno.filt <- metadata[metadata$sampleID %in% colnames(counts.raw),]
pheno.filt <- pheno.filt[match(colnames(counts.raw),pheno.filt$sampleID),]
stopifnot(all(colnames(counts.raw) == pheno.filt$sampleID))

## filter peaks with low counts ##
counts.dge <- DGEList(counts=counts.raw, samples=pheno.filt)
keep <- filterByExpr(counts.dge, group = pheno.filt$fraction)
counts.dge.filt <- counts.dge[keep,]

dim(counts.dge.filt)
dds.filt <- DESeqDataSetFromMatrix(countData = counts.dge.filt,
                                   colData = pheno.filt,
                                   design = ~ fraction)

## normalise counts in peaks with VST transform ##
vst <- vst(dds.filt)
counts.norm <- assay(vst)

## load differential peaks ##
load(paste0(metaDir,"/cellType_diffPeaks_edgeR_5e-2_passS12.narrowPeak.rdata"))
res_sorted <- results[order(results$FDR, decreasing = FALSE),]
most_significant_peaks <- rownames(head(res_sorted, 1000))
counts.dif <- counts.norm[most_significant_peaks,]
counts_to_plot <- counts.dif

## calculate sillhouette coefficient to exclude failed samples
pca_transform <- function(data){
  pca <- prcomp(t(data), center = TRUE, scale. = FALSE)
  
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  plot(var_explained[1:20], type = "b", pch = 19,
       xlab = "PC", ylab = "Proportion of variance explained")
  cum_var <- cumsum(var_explained)
  pc_scores <- pca$x[, 1:5]
  
  return(pc_scores)
}

## Calculate mean silhouette coefficient ##
calc.sil <- function(data, group_labels) {
  
  # Check if data is a matrix or data frame and convert if necessary
  if (!is.matrix(data)) data <- as.matrix(data)
  d <- dist(data)  # Euclidean distance in PC space
  sil <- silhouette(as.integer(as.factor(group_labels)), d)
  
  results <- data.frame(
    sample = rownames(data),
    assigned_cluster = group_labels,
    cluster_number = sil[, "cluster"],
    neighbor_cluster = sil[, "neighbor"],
    sil_width = sil[, "sil_width"]
  )
  
  return(results)
}

counts_to_plot.s <- pca_transform(counts_to_plot)
counts.silh.norm<- calc.sil(counts_to_plot.s,pheno.filt$fraction)
counts.silh.norm <- counts.silh.norm[match(pheno.filt$sampleID,counts.silh.norm$sample),]
pheno.filt$QCS3 <- ifelse(counts.silh.norm$sil_width > silScore, "TRUE", "FALSE")  

# ########################## Fig 1.b ###########################################
# Dotplot of cell-type clustering of samples
# PC1 (x-axis) vs PC2 (y-axis)
# samples coloured by assigned cell-type (color)
# samples marked on pass/fail Stage 3 QC (shape)

pca.df <- prcomp(t(counts_to_plot))
tpm.scores = pca.df$x
tpm.pca<-pca.df$sdev^2/sum(pca.df$sdev^2)
tpm.scores<-tpm.scores[,which(tpm.pca > 0.01)]
percentVar = round(100 * (tpm.pca) )
tpm.scores.df <- data.frame(pc1 = tpm.scores[,1],pc2 = tpm.scores[,2],pc3 = tpm.scores[,3],
                            pc4 = tpm.scores[,4],
                            pc5 = tpm.scores[,5],
                            fraction = as.factor(pheno.filt$fraction),
                            cohort = as.factor(pheno.filt$cohort),
                            sampleCode =as.factor(pheno.filt$sampleCode), 
                            sampleID =as.factor(pheno.filt$sampleID), 
                            gender =as.factor(pheno.filt$gender), passed = pheno.filt$QCS3)

p1b <- ggplot(tpm.scores.df, aes(x=pc1, y=pc2, color =  fraction, shape = passed, label = sampleCode))+
    geom_point(size=2, alpha = 0.7)+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_color_manual(values = col_vals) +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance"), 
         color = "Cell-type", 
         title    = "B", 
         subtitle = "Outliers defined by negative silhouette coefficient", 
         shape = "Pass QC Stage 3") +
    theme_bw(base_size = 11) +
  theme(
      legend.position = "bottom",legend.key.size = unit(1.3, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),
      axis.text = element_text(size =12), plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),ncol = 1),
              colour = guide_legend(override.aes = list(size=6),ncol = 2))

# ########################## Fig 1.c-g ###########################################
# Violin plots of variance explained by each Principal Component (y-axis)
# for each cell-type (x-axis)
# coloured by cell-type (fill)
# marked samples by pass/fail stage 3 QC (shape)

counts_to_plot.m <- melt(tpm.scores.df)
counts_to_plot.m.pc1 <- counts_to_plot.m[counts_to_plot.m$variable == "pc1",]
p1c <- ggplot(counts_to_plot.m.pc1, aes(x=fraction, y = value))+
  geom_jitter(aes(shape = passed))+
  geom_violin(aes(fill = fraction))+
  stat_summary(fun = "median", geom = "point",color = "white",
               position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_fill_manual(values = col_vals) +
    labs(x="",
         y="PC1", 
         fill = "Cell type",  title = "C",
         shape = "Pass Stage 3 QC") +
    theme_bw()+
    theme(
      legend.position = "bottom",legend.key.size = unit(0.6, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),axis.text.x = element_blank(),
      axis.text = element_text(size =12),plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),ncol = 1),
              fill = guide_legend(override.aes = list(size=3),nrow = 1))

counts_to_plot.m.pc2 <- counts_to_plot.m[counts_to_plot.m$variable == "pc2",]
p1d <- ggplot(counts_to_plot.m.pc2, aes(x=fraction, y = value))+
  geom_jitter(aes(shape = passed))+
  geom_violin(aes(fill = fraction))+
  stat_summary(fun = "median", geom = "point",color = "white",
               position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_fill_manual(values = col_vals) +
    labs(x="",
         y="PC2", 
         fill = "Cell type", title = "D",
         shape = "Pass Stage 3 QC") +
    theme_bw()+
    theme(
      legend.position = "bottom",legend.key.size = unit(0.6, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),axis.text.x = element_blank(),
      axis.text = element_text(size =12),plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),ncol = 1),
              fill = guide_legend(override.aes = list(size=3),nrow = 1))

counts_to_plot.m.pc3 <- counts_to_plot.m[counts_to_plot.m$variable == "pc3",]
p1e <- ggplot(counts_to_plot.m.pc3, aes(x=fraction, y = value))+
  geom_jitter(aes(shape = passed))+
  geom_violin(aes(fill = fraction))+
  stat_summary(fun = "median", geom = "point",color = "white",
               position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_fill_manual(values = col_vals) +
    labs(x="",
         y="PC3", 
         fill = "Cell type",  title = "E",
         shape = "Pass Stage 3 QC") +
    theme_bw()+
    theme(
      legend.position = "bottom",legend.key.size = unit(0.6, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),axis.text.x = element_blank(),
      axis.text = element_text(size =12),plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),nrow = 1),
              fill = guide_legend(override.aes = list(size=3),nrow = 1))

counts_to_plot.m.pc4 <- counts_to_plot.m[counts_to_plot.m$variable == "pc4",]
p1f <- ggplot(counts_to_plot.m.pc4, aes(x=fraction, y = value))+
  geom_jitter(aes(shape = passed))+
  geom_violin(aes(fill = fraction))+
  stat_summary(fun = "median", geom = "point",color = "white",
               position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_fill_manual(values = col_vals) +
    labs(x="",
         y="PC4", 
         fill = "Cell type",  title = "F",
         shape = "Pass Stage 3 QC") +
    theme_bw()+
    theme(
      legend.position = "bottom",legend.key.size = unit(0.6, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),axis.text.x = element_blank(),
      axis.text = element_text(size =12),plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),nrow = 1),
              fill = guide_legend(override.aes = list(size=3),nrow = 1))

counts_to_plot.m.pc5 <- counts_to_plot.m[counts_to_plot.m$variable == "pc5",]
p1g <- ggplot(counts_to_plot.m.pc5, aes(x=fraction, y = value))+
  geom_jitter(aes(shape = passed))+
  geom_violin(aes(fill = fraction))+
  stat_summary(fun = "median", geom = "point",color = "white",
               position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_fill_manual(values = col_vals) +
    labs(x="",
         y="PC5", 
         fill = "Cell type",  title = "G",
         shape = "Pass Stage 3 QC") +
    theme_bw()+
    theme(
      legend.position = "bottom",legend.key.size = unit(0.6, 'cm'),
      legend.text = element_text(size =12),legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 14),axis.text.x = element_blank(),
      axis.text = element_text(size =12),plot.title = element_text(size = 16, face = "bold")
    )+ guides(shape = guide_legend(override.aes = list(size=6),nrow = 1),
              fill = guide_legend(override.aes = list(size=3),nrow = 1))

## gather plots ##
p1b_clean <- p1b + theme(legend.position = "none")
p1c_clean <- p1c + theme(legend.position = "none")
p1d_clean <- p1d + theme(legend.position = "none")
p1e_clean <- p1e + theme(legend.position = "none")
p1f_clean <- p1f  + theme(legend.position = "bottom")
p1g_clean <- p1g + theme(legend.position = "none")

figure1 <-
  (p1b_clean |  p1c_clean | p1d_clean) /
  (p1e_clean |  p1f_clean | p1g_clean)
  plot_layout(
    heights = c(1,1),
    nrow = 2
  )

## Save figure to be merged with Figure 1a ##
saveRDS(figure1, "fig1bg_ctCheck.rds")


