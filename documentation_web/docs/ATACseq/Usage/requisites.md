---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# Requisites

To ensure PIPPA can run as expected, please check:

1. The full path to the project folder must be specified as the first argument on the command line for all scripts.
2. The pipeline requires config files (file paths should be specified in config.txt): More info in "Configuration" page
  - config.txt: with variables for bash scripts.
  - config.r with parameters for R scripts
  - condaEnv.yml : file with packages and their version to be installed by conda as an environment.
3. Requires a samples.txt file in the META_DIR (0_metadata) folder with the names of the samples that will be used to run the pipeline. 
4. Requires a sampleSheet.csv file in the META_DIR (0_metadata) folder with information about the samples. Name given/path should be specified as `sampleSheet` in the config.r file. More information below.
5. Samples need to be in the RAWDATADIR folder in the project folder (1_raw)
6. An alphabetically ordered array of cell types to which samples belong to needs to be specified in the *config.txt* file as `CELLTYPES` in order to perform peak calling by cell group (STEP 8). This also required in the *config.r* file. 
7. A Reference folder is needed with the following resources:
  - Reference genome in .fa format (hg38): [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
  - Index for Bowtie2 alignment: [Bowtie index](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)
  - ENCODE blacklist regions: [hg38](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz)
  - Pseudoautosomal regions (PAR) (hg38) [UCSC](https://hgdownload.soe.ucsc.edu/gbdb/hg38/html/description.html)
  - X chromosome gene list (FIRRE and XIST) : chrX	73820656	73852723	XIST	10	- , chrX  131688779	131830604	FIRRE	10	-
  - RefSeq file: [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/)

## Sample sheet

A sample sheet with the metadata about samples is essential to run the ATAC-seq pipeline. This should be located in the META_DIR (0_metadata) folder. The name of this file should be specified in the config files as `SAMPLE_SHEET` and `sampleSheet` in the *config.txt* and *config.r* files, respectively.
Samples must be in the same order as in the samples.txt file. An example of information about a sample is shown below. 

| sequencingBatch | sampleID | sampleCode | cohort | fraction | individualID | age | gender | clinical | vcfID |
| --------------- | -------- | ---------- |------- | -------- | ------------ | --- | ------ | -------- | ----- |
| X	| X_A61_NEUN_S1	| A61_NEUN | Y	| NEUN | 859 | 85	| F |	5	| XXXX_YYYYY_A61 |

- `sequencingBatch`: sequencing batch number
- `sampleID`: name of sample in samples.txt and name of fastq file corresponding to the sample. This is usually formed by the sequencingBatchNumber_sampleID_fraction_SNumber.
- `sampleCode`: ID of sample.
- `cohort`: brain bank or cohort.
- `fraction`: cell type
- `individualID`: ID of the individual to which the sample belongs to. All samples that belong to the same individual must have the same ID.
- `age`
- `gender`: F or M.
- `clinical`: clinical group of samples (optional)
- `vcfID`: ID of imputed array data to perform genotype comparison and confirm sample's identity (optional)

## Other notes

This pipeline has been originally designed to be executed with cell type sorted ATAC-seq datasets. However, with the right modifications, it can also be used on bulk tissue ATAC-seq datasets. Instead of a list of cell types in the `CELLTYPES` variable, use "BULK" for example. This has not been tested yet.