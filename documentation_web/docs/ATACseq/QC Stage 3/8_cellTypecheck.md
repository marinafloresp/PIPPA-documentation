---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 8. Group peak calling

This script performs stage 3 of QC to check cell type of samples by grouping samples per cell type, performing group level peak calling and gets read counts in peaks called. Samples used are those that have passed previous QC stages 1 and 2. Cell fractions to which samples belong to need to be specified in the *config.txt* file.

 `sbatch ./jobSubmission/8_batchRunGroupPeakCalling.sh (project directory) [STEPS] [GROUP]`

##### -scripts executed-
- ./Rscripts/samplesForGroupAnalysis.r (R config file) (GROUP): Creation of a file with a list of samples of cell group input.
- ./subScripts/groupPeaks.sh (GROUP): Peak calling on merged bam file for cell type group or on samples that belong to the cell group.
- ./Rscripts/countsInPeaksGroup.r (R config file): Obtain read counts in peaks called on each cell group and performs differential analysis to find topmost differential peaks.
- ./Rscripts/collateCellTypeCheck.Rmd (R config file): Output results of cell-groups peak calling and counts in peaks to perform cell type check.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `PEAKS` and `MERGE` steps.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `MERGE`: Collates samples bam files and creates a sorted, merged bam file for cell type group.
  - `PEAKS`: Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode. Run only groupPeaks.sh.
  - `COUNTS`: Gets read counts in peaks called at cell group level and analysis for topmost differential peaks. Run only countsInPeaksGroup.r.
  - `CTCHECK` : Collates results from before in a Rmarkdown report and identifies samples that fail cell type identity.