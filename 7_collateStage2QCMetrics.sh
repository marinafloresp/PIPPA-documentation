#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACQCSummaryStage2S7-%A.log
#SBATCH --error=ATACQCSummaryStage2S7-%A.err
#SBATCH --job-name=ATACQCSummaryStage2S7

## ===================================================================================================================##
##                          ATAC-seq pipeline STEP 2: Collate Stage 1 Results                                         ##
## ===================================================================================================================##
## EXECUTION: sbatch ./jobSubmission/7_collateStage2QCMetrics.sh <project-directory> <STEP>                           ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script collates results from the second part of the pipeline: confirmation of correct individual.|| 
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project dir> directory of project, location of config.txt file                                              ||
## $2 -> <STEP> Specify step to run: COLLATE, SUMMARY. Can be combined. Default is to run all                         ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Config.txt file in <project directory>.                                                                          ||
##   SCRIPTS_DIR, PROJECT,FASTQCDIR,CONDA_ENV, CONDA                                                                  ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subscripts                                                            ||
## - Subscripts: progressReportS2.sh                                                                                  ||
## - R subscripts to be in ${RSCRIPTS_DIR} = ./Rscripts                                                               ||
## - R scripts: collateSampleChecks.Rmd                                                                               ||
## - No array job number is required.                                                                                 ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time when script runs
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv ATACQCSummaryStage2S7-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR 
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir: $LOG_DIR

EOF

## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}


## ================ ##
##    VARIABLES     ##
## ================ ##

if [[ $2 == '' ]] 
then
  echo "No step specified. All steps in the script will be run. "
fi

## check step method matches required options and exit if not
if [[ ! $2 == "COLLATE" ]] && [[ ! $2 == "SUMMARY" ]] && [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use COLLATE, SUMMARY or some combination of this as a single string (i.e. MULTIQC,COLLATE)" ; exit 1; }            
fi


## ============ ##
##    STEPS     ##
## ============ ##


## option COLLATE: collate other QC metrics
if [ $# = 1 ] || [[ $2 == 'COLLATE' ]]
then

cat <<EOF

|| Running step 7.1 of ATAC-seq pipeline: COLLATE ||

It will be checked that all samples had gone through all processes and had produced the right outputs.
Output directory is ${META_DIR}

EOF
 
	sh "${SUB_SCRIPTS_DIR}/progressReportS2.sh"
	
fi

## option SUMMARY: Sex check and Genotype check results are collated in Rmarkdown   
if [ $# = 1 ] || [[ $2 == 'SUMMARY' ]]
then

cat <<EOF

|| Running STEP 7.2 of ATAC-seq pipeline: SUMMARY. ||

Results from genotype and sex check will be collated in a Rmarkdown report.
Output directory will be: ${PEAK_DIR}/QCOutput

EOF
 
	Rscript -e "rmarkdown::render('${RSCRIPTS_DIR}/collateSampleChecks.Rmd', output_file='$PEAK_DIR/QCOutput/stage2SummaryStats.html')" "${CONFIGR}"
  
fi

conda deactivate

echo Job finished on:
date -u