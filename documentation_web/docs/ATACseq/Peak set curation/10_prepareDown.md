---
sidebar_position: 3
displayed_sidebar: documentationSidebar

---

# 10. Prepare for downstream analysis

This scripts prepares inputs for downstream analyses, especially counts in peaks matrices. Cell fractions to which samples belong to need to be specified in the *config.txt* file. 

`sbatch ./jobSubmission/9_batchDownstreamAnalyses.sh (project directory) [STEPS] [GROUP] [THRES]`

##### -scripts executed-
- ./Rscripts/curatedPeaksChecks.Rmd (R config file): Output results of curated cell-groups peak calling.
- ./Rscripts/countsInCleanPeaksCT.r (R config file): Obtain read counts in peaks called on each cell group and performs differential analysis to find topmost differential peaks from peak sets used for curation (5e-2). This would be with only samples that passed all QC stages. Note STEP 9.3 should have run before this and created the peak sets.
- ./Rscripts/countsInCuratedPeaks.r (R config file) (GROUP) (THRES): Obtain read counts in peaks called in specified cell type and significance threshold. Note STEP 9.5 PEAKSTR should have run before this and created the peak set.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `COUNTS` step.
- `(THRES)` Select significance threshold of peaks for counts in peaks step 10.3 `COUNTS`. 
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `SUMMARY`: Collates curated peak sets generated across cell types in a Rmarkdown report with  properties of peaks: number, genomic annotation, width, etc. Run curatedPeaksChecks.Rmd only.
  - `CTCLEAN`: Gets read counts in peaks called at cell group level and analysis for topmost differential peaks. Run only countsInCleanPeaksCT.r.
  - `COUNTS`: Gets read counts in peaks called in curated peak set in chosen cell type. Run only countsInCuratedPeaks.r.
