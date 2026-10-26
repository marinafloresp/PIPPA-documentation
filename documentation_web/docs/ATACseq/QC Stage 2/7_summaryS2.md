---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 7. Stage 2 summary

This scripts collates results of Stage 2 of QC pipeline to confirm samples identity. Sex and genotype check results are collated.

`sbatch ./jobSubmission/7_collateStage2QCMetrics.sh (project directory) [STEPS]`

##### -scripts executed-

- ./subScripts/progressReportS2.sh : Identification of samples that have been through stage 2 of QC in order to avoid missing processing samples
- ./Rscripts/collateSampleChecks.Rmd (R config file): Output results from sex and genotype check in a Rmarkdown report.

##### -parameters-
- `--array` : no need to be specified.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `COLLATE`: Collate results from all samples to show progress of each sample through stage 2 of QC pipeline in order to avoid missing steps or samples.
  - `SUMMARY`: Collates stage 2 QC results in a single Rmarkdown report. Run only collateSampleChecks.Rmd.