---
sidebar_position: 2
displayed_sidebar: documentationSidebar

---

# 9. Peak set curation 

This script performs the analysis and curation of peak sets using high-quality samples in order to produce a curated peak set that can be used for downstream analyses. Cell fractions to which samples belong to need to be specified in the *config.txt* file.

 `sbatch ./jobSubmission/9_batchPeakSetCuration.sh (project directory) [STEPS] [GROUP] [THRES]`

##### -scripts executed-

- ./Rscripts/samplesForDownAnalysis.r (R config file) (GROUP): Creation of a file with a list of samples of cell group input.
- ./subScripts/peaksForCuration.sh (GROUP): Peak calling on merged bam file for cell type group or on samples that belong to the cell group.
- ./Rscripts/cutoffAnalysis.r (R config file) (GROUP): Perform segmented analysis to decide optimal peak calling threshold.
- ./subScripts/peaksStrThres.sh (GROUP) (THRES): Peak calling on merged bam file for cell type group or on samples that belong to the cell group with a chosen stringent significance threshold.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `PEAKS`, `PEAKSTR`, and `MERGE` steps.
- `(THRES)` Select significance threshold for peak calling in step 9.5 `STRPEAKS` as suggested by results in `CUTOFF` step 9.4.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `MERGE`: Collates samples bam files and creates a sorted, merged bam file for cell type group.
  - `BW`: Converts merged bam file into bigwig file for genome visualisation.
  - `PEAKS`: Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode. Run only peaksForCuration.sh.
  - `CUTOFF`: Performs segmented analysis to decide an optimal peak calling significance threshold. Run only cutoffAnalysis.r.
  - `PEAKSTR` : Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode with a stringent significance threshold. Run only peaksStrThres.sh.