---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 5. Sex check

This script checks the sex of samples as part of their identity check by creating bam files for each sample containing only reads aligned to sex chromosomes and performing peak calling on these.

  `sbatch --array=<number of batch jobs> ./jobSubmission/5_batchPeaksSexChrs.sh (project directory) [STEPS]`
  
##### -scripts executed-
- ./subScripts/subsetSexChrs.sh (array-number): Subsetting reads of only X and Y chromosomes in bam file for input sample . 
- ./subScripts/sexChrPeaks.sh : Peak calling on the sex chromomes, filtering and read counts using MACS3 in paired-end mode.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. Only for step 5.1 `SPLIT`.
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `SPLIT`: Generate only samples with X and Y chromosomes reads. Run only subsetSexChrs.sh on sample.
  - `PEAKS`: Perform only peak calling on sex chromosomes reads. Run only sexChrPeaks.sh.

##### -requires-
- List of pseudoautosomal regions specified in *config.txt* as `PAR`.
- X chromosome gene list (FIRRE and XIST) specified in *config.txt* as `XCHRBED`.