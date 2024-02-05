### TEMPLATE README.md file (DELETE WHEN COMPLETED) ###
- Replace <text> with details specific to your pipeline
- How to include links: `[text for link](link)`
- How to make bullet points: `* text` or `- text`
- Save README file as .md
###


# < Pipeline title > 


## Pipeline overview

<Include a brief description of:
- the purpose of the pipeline.
- the type of data the pipeline analyses/processes.
- if the pipeline is related to a project(s). >

{e.g. 

**A**ssay for **T**transposase-**A**ccessible **C**hromatin using **seq**uencing (**ATAC-seq**) is a method to investigate the accessibility of choromatin and thus to determine regulatory mechanisms of gene expression. It can help to identify promoter regions and potential enhancers and silencers. ATAC-Seq has become popular for identifying accessible regions of the genome as it’s easier, faster and requires less cells than alternative techniques. 

We have develop a pipeline that process and analyses ATAC-seq data, including pre-processing of samples, alignment, quality control, calling peaks and differential accessibility analysis. The steps to run this pipeline are detailed below.
}


## Requisites of the pipeline

<Include:
- Specify directory to submit the scripts from.
- Location of the config file required.
- Location of the project/data folder: raw data and output. >

{e.g

- All steps should be run from the scripts directory.
- Users require a username folder in the ATACSeq/logFiles directory.
- The name of the project folder must be specified as the first argument on the command line `(project name)`.
- Requires two config files in the ATACSeq/config directory:
  - config.txt with variables for bash scripts
  - config.r with parameters for R scripts
    
}


## Steps

### <Step 1>

< Include:
- Brief description of what is done during this step.
- Each script that is run with a brief description of what it does and their location.
- Basic Command line to execute the script
- Specify the parameters each script takes and whether they are optional (specified as [parameter]) or not (specified as (parameter))
- If optional paremeters, list all the options and their description. >
 
{e.g. 


### 1. Pre-analysis

`sbatch --array="number of batch jobs" ATACSeq/jobSubmission/1_batchRunAlignment.sh (project name) [STEPS]`

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

##### -scripts executed-
- preScripts/fastqc.sh : FastQC for pre-alignment quality control.
- preScripts/fastp.sh : FastP for trimming samples.
- ATACSeq/preprocessing/alignment.sh : Alingment of samples to reference genome using Bowtie2.
- ATACSeq/preprocessing/calcENCODEQCMetrics.sh : ENCODE QC metrics are calculated on aligned samples.
  
##### -parameters-
- `--array`: number of batch jobs, each number matches a sample.
- `(project-name)` the name of the project folder in the main data type data directory
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FASTQC`: Perform only quality control. Run only fastqc.sh on samples.
  - `TRIM`: Perform only trimming on samples. Run only fastp.sh on samples.
  - `ALIGN`: Align samples. Run only alignment.sh on samples.
  - `ENCODE`: Calculation of ENCODE QC metrics.

}

### <Step 2>

### <Step 3>

### <Step ...>


## Contributors (Optional) 

< Include who contributes to the development of the pipeline. >
	
## Publications (Optional)

<Include related publications to the project/pipeline.>

## References/Resources (Optional)

<Include references used to build the pipeline, useful resources or guidance used.>

{e.g. 

The development of this pipeline has been largely based on the approach described in: [Yan, F., Powell, D.R., Curtis, D.J. et al. From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis. Genome Biol 21, 22 (2020). https://doi.org/10.1186/s13059-020-1929-3](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3#citeas)

}
