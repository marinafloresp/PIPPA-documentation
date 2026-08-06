## ===============================================================================================================##
##             ATAC-seq pipeline STEP 9.0: Sample list for downstream analyses                                    ##
## ===============================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/samplesForDownAnalysis.r <project> <cell-group>                ||
## - execute from scripts directory                                                                               ||
##                                                                                                                ||
## DESCRIPTION: This scripts creates a file with samples that passed all stages of the QC and selects             ||
## samples that belong to an input cell type to perform peak calling.                                             ||
##                                                                                                                ||
## INPUTS:                                                                                                        ||
## - <project> : project on which analysis is being run                                                           ||
## - <ell-group>: cell fraction of samples to select for group peak calling                                       ||
##                                                                                                                ||
## OUTPUTS:                                                                                                       ||
## - downAnalysisSamples_<cell-group>.txt" in 0_metadata folder                                                   ||
##                                                                                                                ||
## REQUIRES:                                                                                                      ||
## - R/4.2.1-foss-2022a                                                                                           ||
## - A csv file in the metadata folder with the summary of QC stages 1 and 2: passStage3Status.csv. This is       ||
## produced at STEP 8.2 CTCHECK                                                                                   ||
## ===============================================================================================================##

args <- commandArgs(trailingOnly=TRUE)
configFile <-args[1]
source(configFile)
cf <- args[2]

passAllQC <- read.csv(file.path(allQCresults), stringsAsFactors = FALSE, strip.white = TRUE)
samples <- passAllQC[passAllQC$fraction == cf & passAllQC$PASSALL == TRUE,]$sampleID
samples <- na.omit(samples)
write.table(samples, file = file.path(paste0(metaDir, "/downAnalysisSamples_",cf,".txt")), row.names = FALSE,quote = FALSE, col.names=FALSE)
