#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACmacsSimulation-%A.log
#SBATCH --error=ATACmacsSimulation-%A.err
#SBATCH --job-name=ATACmacsSimulation

## ===================================================================================================================##
##                                      ATAC-seq analysis: MACS analysis                                              ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./jobSubmission/macsSimulation.sh <project directory> <STEP> <GROUP>                            ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This scripts recreates the analysis done to explore MACS' peak calling behaviour.                      ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <project directory> directory of project, location of config.txt file                                        ||
## $2 -> <STEP> Specify step to run: MERGE, SUB, PEAKS, COUNT.                                                        ||
## $3 -> <GROUP> Specify cell group to run previous step(s).                                                          ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${META_DIR}/downAnalysisSamples_${GROUP}.txt that lists sample names.                                    ||
## - Config.txt file in <project directory>.                                                                          ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, ALIGNED_DIR, PEAK_DIR,            ||
##   SCRIPTS_DIR, PROJECT,CONFIGR,RSCRIPTS_DIR,CONDA_ENV, CONDA                                                       ||
## - For modules or references required, please refer to each subscript run in this script.                           ||
## - A conda environment setup with several modules: samtools, MACS3                                                  ||
## - Subscripts to be in ${SUB_SCRIPTS_DIR} = ./subScripts                                                            ||
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
mv ATACmacsSimulation-${SLURM_ARRAY_JOB_ID}* $LOG_DIR

cat <<EOF

Loading config file for project:  ${PROJECT}
Project directory is:  $MAIN_DIR
Script is running from directory:  ${SCRIPTS_DIR}
Log files will be moved to dir:  $LOG_DIR

EOF

## Activate conda environment with packages/modules
source ${CONDA} 
conda activate ${CONDA_ENV}
#conda activate newSam

## ================ ##
##    VARIABLES     ##
## ================ ##

export GROUP=$3
export NUM=$4

if [[ ! ${GROUP} == '' ]]
then
  echo "Chosen cell fraction is ${GROUP}"

  ## Retrieve samples from the same cell fraction that passed QC
  mapfile -t SAMPLES < ${META_DIR}/downAnalysisSamples_${GROUP}.txt
  sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
  echo "Number of samples found in ${GROUP} is:" """${#SAMPLES[@]}"""
else
  echo "No cell group specified."
fi

if [[ ! $2 == '' ]]
then
  echo "Chosen step is $2"
else
  echo "No step specified."
fi

## ========= ##
##   STEPS   ##
## ========= ##

if [ $# == 1 ] || [[ $2 == 'MERGE' ]]
then

  OUPUTDIR="${ALIGNED_DIR}/merged/passAllQC/"
  mkdir -p "${OUPUTDIR}"
  
  cat <<EOF

|| Running STEP 1 of MACS analysis: MERGE ||

A merged bam file will be created with all samples that passed all QC for the chosen cell-type.
Output directory is ${OUPUTDIR}
Aligned files taken from ${ALIGNED_DIR}

EOF

  bamFiles=()
  for sampleName in "${SAMPLES[@]}"; do
    bamFiles+=("${ALIGNED_DIR}/${sampleName}.filt.nodup.bam")
  done
  
  echo "Number of aligned bam files found in ${ALIGNED_DIR} is: ${#bamFiles[@]}"

  samtools merge -f "${OUPUTDIR}/${GROUP}.merged.bam" "${bamFiles[@]}"
  samtools sort -@ 4 -o "${OUPUTDIR}/${GROUP}.sorted.merged.bam" "${OUPUTDIR}/${GROUP}.merged.bam"
  echo "Merged sample is sorted"
  samtools view -H "${OUPUTDIR}/${GROUP}.sorted.merged.bam" > "${OUPUTDIR}/${GROUP}.sorted.merged_header.sam"
  echo "Merged sample with samples that belong to ${GROUP} has been created."
  rm "${OUPUTDIR}/${GROUP}.merged.bam"
  
fi

if [ $# = 1 ] || [[ $2 == 'COUNT' ]]
then

  BAM_FILE="${ALIGNED_DIR}/merged/passAllQC/${GROUP}.sorted.merged.bam"

  cat <<EOF

|| Running STEP 2 of MACS analysis: COUNT ||

Count number of total reads in merged bam file created in STEP 1 MERGE.
Bam file read is ${BAM_FILE}

EOF

  total_reads=$(samtools view -c "$BAM_FILE")
  total_pairs=$(( total_reads / 2 ))
  echo "Total number of reads: $total_reads"
  echo "Total read pairs: $total_pairs"
  
fi

if [ $# = 1 ] || [[ $2 == 'DOWN' ]]
then
  BAM_FILE="${ALIGNED_DIR}/merged/passAllQC/${GROUP}.sorted.merged.bam"
  OUTPUTDIR="${ALIGNED_DIR}/numReadsSubsets/${GROUP}/"
  mkdir -p "${OUTPUTDIR}"
  
  SCALE=10000000
  readPairsToSample=$((${NUM} * SCALE))
  
  cat <<EOF

|| Running STEP 3 of MACS analysis: DOWN ||

Create subset bed file with randomly sampled reads from merge bam file created in STEP 1 MERGE.
Note input number should be read pairs. Selected number is ${NUM}, which will be multiplied by ${SCALE}.
Total number of read pairs downsampled will be ${readPairsToSample}
Bam file read is ${BAM_FILE}
Output directory is ${OUTPUTDIR}

EOF
  
  macs3 randsample -i ${BAM_FILE} -f BAMPE --seed 2 --outdir ${OUTPUTDIR} -o "${OUTPUTDIR}/${GROUP}_downSamp_${NUM}.bed" -n ${readPairsToSample}

  echo "Downsampled bed file with $(wc -l "${OUTPUTDIR}/${GROUP}_downSamp_${NUM}.bed" | cut -f1 -d' ') reads have been created "
  
fi


if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
  
  INPUTDIR="${ALIGNED_DIR}/numReadsSubsets/"
  BEDFILE="${INPUTDIR}/${GROUP}_downSamp_${NUM}.bed"
  MACS_GROUP=${THRES}
  
  PEAK_DIR_GROUPS="${PEAK_DIR_GROUPS}/numSamplesPeaks/"
  PEAKS_NAME="${GROUP}_${MACS_GROUP}_${NUM}down"

  SCALE=10000000
  readPairsToSample=$((${NUM} * SCALE))
  
  mkdir -p "${PEAK_DIR_GROUPS}"
  
 cat <<EOF

|| Running STEP 4 of MACS analysis: PEAKS ||

Peak calling will be performed on downsampled bed file created in STEP 3 DOWN.
Bed file created with ${readPairsToSample} downsampled reads will be used: "${INPUTDIR}/${GROUP}_downSamp_${NUM}.bed"
Output directory is ${PEAK_DIR_GROUPS}
Peak file name is ${PEAKS_NAME}
Q-value threshold is: ${MACS_GROUP}

EOF

  macs3 callpeak -t ${BEDFILE} --outdir ${PEAK_DIR_GROUPS} -n ${PEAKS_NAME} -f BEDPE -g 2.9e9 -q $MACS_GROUP --keep-dup all --nomodel 

  ## exclude peaks aligned to blacklist regions and peaks called in chr X and Y
  bedtools intersect -v -a ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.narrowPeak -b ${BLACKLIST} \
  	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  	| awk '!/^(chrY|chrX)/' > ${PEAK_DIR_GROUPS}/${PEAKS_NAME}.narrowPeak.filt
   
  rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.narrowPeak
  rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.xls 
  
fi
if [ $# = 1 ] || [[ $2 == 'PEAKS' ]]
then
   
  INPUTDIR="${ALIGNED_DIR}/numReadsSubsets/${GROUP}/"
  BEDFILE="${INPUTDIR}/${GROUP}_downSamp_${NUM}.bed"
  MACS_GROUP=$5
  
  PEAK_DIR_GROUPS="${PEAK_DIR_GROUPS}/numSamplesPeaks_test/${GROUP}/"
  PEAKS_NAME="${GROUP}_${MACS_GROUP}_${NUM}down"
  
  SCALE=10000000
  readPairsToSample=$((${NUM} * SCALE))
  
  mkdir -p "${PEAK_DIR_GROUPS}"
  
 cat <<EOF

|| Running STEP 4 of MACS analysis: PEAKS ||

Peak calling will be performed on downsampled bed file created in STEP 3 DOWN.
Bed file created with ${readPairsToSample} downsampled reads will be used: "${INPUTDIR}/${GROUP}_downSamp_${NUM}.bed"
Output directory is ${PEAK_DIR_GROUPS}
Peak file name is ${PEAKS_NAME}
Q-value threshold is: ${MACS_GROUP}

EOF

  macs3 callpeak -t ${BEDFILE} --outdir ${PEAK_DIR_GROUPS} -n ${PEAKS_NAME} -f BEDPE -g 2.9e9 -q $MACS_GROUP --keep-dup all --nomodel 

  ## exclude peaks aligned to blacklist regions and peaks called in chr X and Y
  bedtools intersect -v -a ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.narrowPeak -b ${BLACKLIST} \
  	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  	| awk '!/^(chrY|chrX)/' > ${PEAK_DIR_GROUPS}/${PEAKS_NAME}.narrowPeak.filt
   
  rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.narrowPeak
  rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME}_peaks.xls 
  
fi

conda deactivate
 
echo Job finished on:
date -u

