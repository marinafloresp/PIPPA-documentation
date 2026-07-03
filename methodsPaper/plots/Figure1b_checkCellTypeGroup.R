## Figure 1.b: Stage 3 QC Cell-type check

library(edgeR)
library(ggplot2)
library(DESeq2)
library(cluster)

## Load metadata ##
metadata <- read.table(paste0(metaDir,"/methodsPaper/sampleSheet_50bp.csv"), header = TRUE, sep = ',', stringsAsFactors = FALSE)
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
dim(counts_to_plot)

## Plot counts in peaks (before filtering) ##
col_vals <- c("Sox10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")
## Plot PCA with first two principal components
pca_single <- function(counts,pheno_sub, col){
  pca <- prcomp(t(counts))
  tpm.scores = pca$x
  tpm.pca<-pca$sdev^2/sum(pca$sdev^2)
  tpm.scores<-tpm.scores[,which(tpm.pca > 0.01)]
  percentVar = round(100 * (tpm.pca) )
  tpm.scores.df = cbind(pc1 = tpm.scores[,1],pc2 = tpm.scores[,2],
                        pc3 = tpm.scores[,3],pc4 = tpm.scores[,4]
                        ,pheno_sub)
  plot1 <- ggplot(tpm.scores.df, aes(x=pc1, y=pc2, color = fraction,label= sampleCode))+
    geom_point(size=4, alpha = 0.7)+
    labs(x=paste0("PC1 (", percentVar[1], "%)"),y=paste0("PC2 (", percentVar[2], "%)"),  
         title = "", color = "Cell type") +
    scale_color_manual(values = col_vals) +
    theme_bw() +
    theme(
      legend.position = "bottom",legend.key.size = unit(1.3, 'cm'),
      legend.text = element_text(size =20),legend.title = element_text(size = 20),
      axis.title = element_text(margin = margin(r = 15), size = 20),
      axis.text = element_text(size =20), 
    )+ guides(shape = guide_legend(nrow = 2))
  
  return(plot1)
}
p1 <- pca_single(counts_to_plot, pheno.filt)
p1


## calculate sillhouette coefficient to exclude failed samples
pca_transform <- function(data){
  pca <- prcomp(t(data), center = TRUE, scale. = FALSE)
  
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  plot(var_explained[1:20], type = "b", pch = 19,
       xlab = "PC", ylab = "Proportion of variance explained")
  cum_var <- cumsum(var_explained)
  
  # Number of PCs needed to reach 70%
  n_pcs_70 <- which(cum_var >= 0.70)[1]
  pc_scores <- pca$x[, 1:n_pcs_70]
  
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

## plot samples with assigned sillhouette score
pca_qcs3 <- function(counts,pheno_sub){
  pca <- prcomp(t(counts))
  tpm.scores = pca$x
  tpm.pca<-pca$sdev^2/sum(pca$sdev^2)
  tpm.scores<-tpm.scores[,which(tpm.pca > 0.01)]
  percentVar = round(100 * (tpm.pca) )
  tpm.scores.df <- data.frame(pc1 = tpm.scores[,1],pc2 = tpm.scores[,2],fraction = as.factor(pheno_sub$fraction),
                              cohort = as.factor(pheno_sub$cohort),sampleCode =as.factor(pheno_sub$sampleCode), 
                              sampleID =as.factor(pheno_sub$sampleID), gender =as.factor(pheno_sub$"gender"), passed = pheno_sub$QCS3)
  
  plot1 <- ggplot(tpm.scores.df, aes(x=pc1, y=pc2, color =  fraction,
                                     shape = passed, label = sampleID))+
    geom_point(size=4, alpha = 0.7)+
   # stat_ellipse(geom = "polygon", alpha = 0.15, color = NA, level = 0.95)+
  scale_shape_manual(values = c("FALSE" = 3, "TRUE" = 16)) + 
    scale_color_manual(values = col_vals) +
   # scale_fill_manual(values = col_vals) +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance"), 
         color = "CellType", 
         title    = "Cell Type Check: PCA of top 1000 cell differential ATAC peaks",
         subtitle = "Outliers defined by negative silhouette coefficient", 
         shape = "Pass QC Stage 3") +
    
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1.5, 'cm'),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 12),
      axis.title = element_text(margin = margin(r = 15), size = 15),
      axis.text = element_text(size =15), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size =15)
    )+ guides(shape = guide_legend(override.aes = list(size=6),ncol = 1),
              colour = guide_legend(override.aes = list(size=6),ncol = 2))
  
  
  return(plot1)
}

counts.silh.norm <- counts.silh.norm[match(pheno.filt$sampleID,counts.silh.norm$sample),]
pheno.filt$QCS3 <- ifelse(counts.silh.norm$sil_width < 0, "FALSE", "TRUE")
table(pheno.filt$QCS3)
p2 <- pca_qcs3(counts_to_plot,pheno.filt)

## save new samplesheet with passed Stage 3 QC samples
stage3QC <- pheno.filt
qc2 <- read.csv(file.path(paste0(metaDir,"/methodsPaper/passStage2Status_newMetrics.csv")), sep = ",",stringsAsFactors = FALSE)
qc2$QCPASS <- as.factor(qc2$QCPASS)
stage3QC$QCS2 <- qc2[match(stage3QC$sampleID,qc2$sampleID),]$QCPASS
qc1 <- read.csv(file.path(paste0(metaDir,"/methodsPaper/passStage1Status_newMetrics.csv")), sep = ",",stringsAsFactors = FALSE)
qc1$QCPASS <- as.factor(qc1$QCPASS)
stage3QC$QCS1 <- qc1[match(stage3QC$sampleID,qc1$sampleID),]$QCPASS

stage3QC$PASSALL <- ifelse(stage3QC$QCS1 == TRUE & stage3QC$QCS2 == TRUE & stage3QC$QCS3 == TRUE, TRUE, FALSE)
table(stage3QC$PASSALL)

write.csv(stage3QC, file = "passStage3Status_newMetrics.csv", row.names = FALSE)

##list of samples
samplesPassAll <- stage3QC[stage3QC$PASSALL == TRUE,]$sampleID

write.table(samplesPassAll,file = "passStage3SampleList.txt", row.names = FALSE, col.names = FALSE)


