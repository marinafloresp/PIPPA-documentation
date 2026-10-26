---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 3. Samples peak calling

This script performs the core part of the ATACseq pipeline: calling peaks. In this stage, peak calling is performed at sample level using the paired-end mode of MACS3. After this, FRIP statistics are calculated and collated in a single csv file.

  `sbatch --array=<number of batch jobs> ./jobSubmission/3_batchRunPeakCalling.sh (project directory) [STEPS]`

##### -scripts executed-

- ./subScripts/samplePeaks.sh (array-number): MACS version 3 peak calling with BAM files paired-end reads
    * it then filters the peaks to exclude those that overlap with blacklisted regions and sex chromosomes.
    * peaks are sorted by chromosome.
- ./subScripts/collateCalcFrip.sh (array-number): Calculation of fraction of reads in peaks, number of reads and peaks for peak calling at sample level
    * a single --array job number should be specified (any number).

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample. No needed for `FRIP`.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `PEAKS`: Perform only peak calling on samples. Run only samplePeaks.sh on sample.
  - `FRIP`: Calculate peak calling results for all samples. Run only collateCalcFrip.sh on all samples found.