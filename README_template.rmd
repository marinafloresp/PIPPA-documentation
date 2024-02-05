### TEMPLATE README.md file (DELETE WHEN COMPLETED) ###
## Replace <text> with details specific to your pipeline
## How to include links: [text for link](link)
## How to make bullet points: * text
## Save README file as .md

# <Pipeline title>

## Pipeline overview

<Include a brief description of:
- the purpose of the pipeline.
- the type of data the pipeline analyses/processes.
- if the pipeline is related to a project(s). >

## Requisites of the pipeline

<Include:
- Specify directory to submit the scripts from.
- Location of the config file required.
- Location of the project/data folder: raw data and output. >

## Steps

### <Step 1>

< Include:
- Brief description of what is done during this step.
- Each script that is run with a brief description of what it does and their location.
- Basic Command line to execute the script
- Specify the parameters each script takes and whether they are optional (specified as [parameter]) or not (specified as (parameter))
- If optional paremeters, list all the options and their description. >
 
{e.g. 1. Data pre-processing

* 1. `sbatch --array="number of batch jobs" ATACSeq/jobSubmission/1_batchRunAlignment.sh (project name) [STEPS]`
Performs first pre-analysis: pre-alignment QC, alignment and post-alignment QC.
	- executes fastqc.sh
	- executes fastp.sh
	- executes ATACSeq/preprocessing/1_alignment.sh
	- executes ATACSeq/preprocessing/2_calcENCODEQCMetrics.sh ${sampleID}_sorted_chr1.bam
##### -parameters-
	- `--array`: number of batch jobs, each number matches a sample.
	- (project-name) the name of the project folder in the main data type data directory
	-optional-
	- [STEPS] option of either FASTQC, TRIM, ALIGN or ENCODE. The first runs only fastqc on the samples, the second trims, the third runs only the alignment and the last calculates the ENCODE QC metrics. 
	Option to combine steps with desired steps included as single string, i.e. FASTQC,TRIM. Default if left blank is to run all steps.
}

### <Step 2>

### <Step 3>

### <Step ...>


## Contributors (Optional)

<Include people that contributes to development of the pipeline>

## Publications (Optional)

<Include related publications to the project/pipeline.>

## References/Resources (Optional)

<Include references used to build the pipeline, useful resources or guidance used.>
