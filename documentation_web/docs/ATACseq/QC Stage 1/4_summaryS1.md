---
sidebar_position: 4
displayed_sidebar: documentationSidebar

---

# 4. Stage 1 summary

This scripts collates results of Stage 1 of QC pipeline to confirm data quality. It uses MultiQC to collate the output of FastQC, Bowtie2 alignment, ATAC-seq QC metrics and sample peak calling results. 

  `sbatch ./jobSubmission/4_collateStage1QCMetrics.sh (project directory) [STEPS]`

##### -scripts executed-
- ./subScripts/progressReport.sh : Identification of samples that have been successful at each stage of the processing pipeline and for each fastq file, how far through the process it has progressed.
- ./subScripts/countMTcollateFS.sh : Counting of the number of reads aligned to MT chromosome and collates flagstat summary of- aligned sorted reads.
- ./Rscripts/collateDataQualityStats.Rmd (R config file): Output Rmarkdown report summarising stage 1 QC metrics: raw reads metrics, trimming, alignment and peak calling.

##### -parameters-
- `--array` : no need to be specified.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `MULTIQC`: Collates fastQC, trimming and alignment statistics in a single report. Run multiQC to collate fastQC, trimming and alignment results.
  - `COLLATE`: Collate results from all samples to show progress of each sample through the pipeline in order to avoid missing steps. Run only progressReport.sh and countMTcollateFS.sh.
  - `SUMMARY`: Collates stage 1 QC and peak calling results in a single Rmarkdown report. Run only collateDataQualityStats.Rmd.