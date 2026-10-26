---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# 6. Genotype check

This script checks the genotype of the samples as part of their identity check by detecting possible DNA contamination in order to ensure high quality sequence reads. If any contamination is detected, possible swaps will be suggested. 

 `sbatch --array=<number of batch jobs> ./jobSubmission/6_batchRunGenotypeCheck.sh (project directory) [STEPS]`

##### -scripts executed-
- ./subScripts/compareBamWithGenotypes.sh (array-number): Preparation of bam file for comparison of assigned genotype of each sample and runs verifyBamID 
  * requires a file in METADATA folder called matchedVCFIDs.txt which lists the samples with their matched vcfID.
- ./subScripts/searchBestGenoMatch.sh : Output a summary of stats from previous step and finds any sample that might be contaminated.
  * Contaminated samples will go to a created file potentialSwitches.txt. If this file is not empty, searchBestGenoMatch.sh is executed and an alternative Genotype search is done for the sample.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. 
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. GENCHECK,COMPARE. Default if left blank is to run all of them.
  - `COMPARE`: Compare only samples genotype with matched VCF file genotype. Run only compareBamWithGenotypes.sh on sample.
  - `SWITCH`: Indentify only mismatched genotype samples and output potential switches and seach for alternative genotype. Run only searchBestGenoMatch.sh.