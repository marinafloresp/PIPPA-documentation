---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 1. Pre-analysis

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

 `sbatch --array=<number of batch jobs> ./jobSubmission/1_batchRunPreAnalysis.sh (project directory) [STEPS]`

##### -scripts executed-
- ./subScripts/fastQC.sh (array-number): FastQC for pre-alignment quality control.
- ./subScripts/fastp.sh (array-number): FastP for trimming samples.
- ./subScripts/alignment.sh (array-number): Alingment of samples to reference genome using Bowtie2.
- ./subScripts/calcENCODEQCMetrics.sh (array-number): ENCODE QC metrics are calculated on aligned samples.
  
##### -parameters-
- `--array`: number of batch jobs, each number matches a sample, starting from 0.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FASTQC`: Perform only quality control. Run only fastQC.sh on sample.
  - `TRIM`: Perform only trimming on samples. Run only fastp.sh on sample.
  - `ALIGN`: Align samples. Run only alignment.sh on sample.
  - `ENCODE`: Calculation of ENCODE QC metrics. Run only calcENCODEQCMetrics.sh on sample.

##### -requires-
- file with adapters to be trimmed at `STEP 1.2 TRIM`. Should be in fasta format .fa and specified in *config.txt* file as `ADAPTERS_FILE`. 