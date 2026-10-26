---
sidebar_position: 1
displayed_sidebar: documentationSidebar
---

# 2. ATAC-seq QC metrics

Performs analysis of ATAC-seq-specific characteristics, including fragment size distribution and TSS enrichment scores.

`sbatch --array=<number of batch jobs> ./jobSubmission/2_batchCalcQCMetrics.sh  (project directory) [STEPS] (GROUP_SAMPS)`

##### -scripts executed-  
- ./Rscripts/fragmentDistribution.r (R config file) (array-number) : Calculation of fragment distribution of samples specified by array. 
- ./subScripts/prepareBigWig.sh (array-number) : Coversion of sample bam file to bigWig format and shift.
- ./subScripts/calcTssEnrich.sh (array-number) : Calculation of TSS enrichment score of sample.
- ./subScripts/collateTssScores.sh (array-number) : Collation of TSS enrichment score of all samples.

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample, starting from 0. Only required for BW and SCORE.
- `(project-name)` project's directory.
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FRAG`: generate the fragment distribution and calculate some summary statistics. Run fragmentDistribution.r on group of 10 samples.
  - `BW`: Transform bam file to bigwig format and shift. Run only prepareBigWig.sh on sample.
  - `SCORE`: Calculate TSS enrichment score. Run only calcTssEnrich.sh on sample.
  - `COLLATE`: Collate TSS enrichment scores of all samples. Run only collateTSSscores.
- `GROUP_SAMPS` : for step `FRAG` , should be the number of samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples each.

##### -requires-
- RefSeq file for TSS scores calculation. Should be specified as `REFTSS` in the *config.txt* file. 