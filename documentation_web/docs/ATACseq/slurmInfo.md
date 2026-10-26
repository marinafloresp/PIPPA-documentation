---
sidebar_position: 6 
---

# SLURM information

**PIPPA** has been developed to be executed using SLURM workload manager, 
so main steps scripts would need to be run through SLURM's `sbatch` command. Examples are given at each step page.

## `#SBATCH` configuration

Every main step scripts (those in _jobSubmission_ folder, starting with a number), starts with a section with several comment lines that start with `#SBATCH`. These lines include the configuration read by SLURM in order to run the script. 

Example below: 
```
#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACAlignment-%A_%a.o
#SBATCH --error=ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment
```
You can change the configuration in each main script or adding these variables when running the script with the `sbatch` command: 
` sbatch --array=0 --job-name=ATACAlignment_test ./jobSubmission/1_batchRunPreAnalysis.sh ./project1 FASTQC`

## Log files

Whenever a main step script is run, it will generate error and output log files. If the step is run globally (for all samples), single error and output log files will be generated (1 .err and 1 .log files). If the step is run in parallel for each sample, like sample-peak calling, then two files will be generated for job submitted (one per sample). Desired location to get these log files should be specified as the `LOG_DIR` variable in the _config.txt_ file.

## Arrays

Some steps of the pipeline have been designed to be run on parallel for each sample, so one job would be submitted per sample. These would be run as a SLURM array (parallelisation). For the steps that require these, scripts need to be run through `sbatch --array=x-y[script] `.

For example:

` sbatch --array=0-5 ./jobSubmission/1_batchRunPreAnalysis.sh ./project1 FASTQC`

In that case, FASTQC output would be generated for samples 2 to 6, as SLURM array starts from 0. 

When the step does not need a SLURM array as it a summary/collating step, you can ommit the `--array` parameter or use `--array=0`.
The requirements of each step has been clearly documented in their corresponding sections.