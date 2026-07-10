#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACcalcQCS2-%A_%a.log
#SBATCH --error=ATACcalcQCS2-%A_%a.err
#SBATCH --job-name=ATACcalcQCS2

## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 2: Post alignment processing                                    ##
## ===================================================================================================================##
## EXECUTION: sbatch --array=<sample-index> ./jobSubmission/2_batchCalcQCMetrics.sh <project directory> <STEP>        ||
## - execute from pipeline's main scripts directory                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script performs the second step in the ATAC-seq pipeline: ATAC-seq-specific metrics calculation. ||
## Includes the classification of aligned reads into fragment sizes and related statistics (FRAG).Then pre-processing ||
## of aligned reads are done (BW) to calculate TSS enrichment score (SCORE) and collate results (COLLATE)             ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## --array -> Number of jobs to run. Will select sample(s) corresponding to the number(s) input. Needed for steps BW, ||
##  and SCORE.                                                                                                        || 
## - $1 -> <project directory> : directory to config file for project                                                 ||
## - $2 -> <STEP> Specify step to run: FRAG, BW, SCORE or COLLATE. Needs to be specified.                             ||
## - $3 -> <GROUP_SAMPS> Only to be specified for FRAG step. It splits the samples into groups of 10 to run in        ||
##  parallel, so batch number should be the number of samples divided by 10. e.g. If there are 50 samples, batch      ||
##  number should be 0-5, producing 5 batches of 10 samples each. If the number of samples is less than 10, set batch ||
##  number to 0.                                                                                                      ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/samples.txt that lists sample names.                                                         ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: MAIN_DIR, SCRIPTS_DIR, LOG_DIR, ALIGNED_DIR, RSCRIPTS_DIR,     ||
##   CONDA, CONDA_ENV                                                                                                 ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, R, samstats,bedtools, bowtie2 and fastp                ||
## - R version should be > 4.2 in conda environment.                                                                  ||
## - fragmentDistribution.r file in ${RSCRIPTS_DIR} = ./Rscripts directory                                            ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
## - Subscripts: prepareBigWig.sh, calcTssEnrich.sh, collateTSSScores.sh                                              ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

source "${1}/config.txt" || { echo "No project directory specified or could not be found." ; exit 1; }

## Log files directory
LOG_DIR=${LOG_DIR}/${USER}/${SLURM_ARRAY_JOB_ID}
mkdir -p $LOG_DIR
mv ATACcalcQCS2-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir:  $LOG_DIR

EOF


## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}

## ================ ##
##    VARIABLES     ##
## ================ ##

if [[ $2 == '' ]] 
then
  { echo "No step specified. All steps in the script will be run. "; exit 1; }
fi

## check step method matches required options and exit if not
if [[ ! $2 == "FRAG" ]] && [[ ! $2 == "BW" ]] && [[ ! $2 == "SCORE" ]] && [[ ! $2 == "COLLATE" ]] &&  [[ ! $2 == '' ]];
then 
    { echo "Unknown step specified. Please use FRAG, BW, SCORE, COLLATE" ; exit 1; }            
fi

## check array specified and exit if not for BW and SCORE steps
if [[ ( $2 == "BW" || $2 == "SCORE" ) && ! -v SLURM_ARRAY_TASK_ID ]];
then 
    { echo "Invalid -array number specified for $2 STEP. Each number corresponds to a sample from the sample list provided." ; exit 1; }
fi

## check group of samples specified and exit if not for FRAG step
if [[ $2 == "FRAG" && $3 == '' ]]
then
  { echo "No group of samples specified for STEP 2.1 FRAG. "; exit 1; }
fi

## Collate all samples
mapfile -t SAMPLESIDS < ${META_DIR}/samples.txt

echo "Total number of samples: "" ""${#SAMPLESIDS[@]}"""

##  ==================  ##
##     CALCQCMetrics    ##
##  ==================  ##

if [ $# == 1 ] || [[ $2 == 'FRAG' ]]
then

  mkdir -p ${ALIGNED_DIR}/QCOutput
  
  GROUP_SAMPS=$3

cat <<EOF

|| Running STEP 2.1 of ATAC-seq pipeline: Post alignment processing (QC metrics and fragment distribution). ||

Calculating QC metrics and fragment distribution for samples in group of samples ${GROUP_SAMPS} 
Output directory is ${ALIGNED_DIR}/QCOutput/

EOF
  
  Rscript "${RSCRIPTS_DIR}/fragmentDistribution.r" ${CONFIGR} ${GROUP_SAMPS} 
  
  if [[ ! -f "${ALIGNED_DIR}/QCOutput/FSD_${GROUP_SAMPS}.pdf" ]] && [[ ! -f "${ALIGNED_DIR}/QCOutput/FragmentDistribution_${GROUP_SAMPS}.rdata" ]]
  then 
    { echo "Fragment distribution and other QC metrics of samples in batch ${GROUP_SAMPS} could not be calculated. Please make sure STEP 1.3 ALIGN was successfully run."  ; exit 1;}
  else
    echo "Fragment distribution and other QC metrics of samples in batch ${GROUP_SAMPS} were successfully calculated."
  fi
fi

if [ $# == 1 ] || [[ $2 == 'BW' ]]
then

  export TSS_DIR="${ALIGNED_DIR}/bigWig/tssEnrich/"
  mkdir -p "${TSS_DIR}"
  
  ##Select sample using --array number
  sampleName=${SAMPLESIDS[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample selected is ${sampleName}"

cat <<EOF

|| Running STEP 2.2 of ATAC-seq pipeline: Prepare file for TSS calculation. ||

Transforming sample ${sampleName} to bigWig format to calculate TSS enrichment score. 
Output directory is ${TSS_DIR}

EOF

  sh "${SUB_SCRIPTS_DIR}/prepareBigWig.sh" ${sampleName}
  
fi


if [ $# == 1 ] || [[ $2 == 'SCORE' ]]
then

  export TSS_DIR="${ALIGNED_DIR}/bigWig/tssEnrich/"
  mkdir -p "${TSS_DIR}"
  
  ##Select sampleName using --array number
  sampleName=${SAMPLESIDS[${SLURM_ARRAY_TASK_ID}]}
  echo "Sample selected is ${sampleName}"

cat <<EOF

|| Running STEP 2.3 of ATAC-seq pipeline: Calculate TSS enrichment. ||

Calculate TSS enrichment score for sample ${sampleName}.
Output directory is ${TSS_DIR}

EOF

  sh "${SUB_SCRIPTS_DIR}/calcTssEnrich.sh" ${sampleName}

fi

if [ $# == 1 ] || [[ $2 == 'COLLATE' ]]
then

  export TSS_DIR="${ALIGNED_DIR}/bigWig/tssEnrich/"
  mkdir -p ${ALIGNED_DIR}/QCOutput  
  
cat <<EOF

|| Running STEP 2.4 of ATAC-seq pipeline: Collate TSS enrichment scores. ||

TSS enrichment scores for all samples are collated in a single file.
Output directory is ${ALIGNED_DIR}/QCOutput 

EOF

  sh "${SUB_SCRIPTS_DIR}/collateTSSscores.sh" 

fi

conda deactivate
 
echo Job finished on:
date -u

