#!/bin/bash
## ===================================================================================================================##
##                    ATAC-seq pipeline STEP 9.2: Cell group Peak Calling (case-control)                              ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/caseControlPeaks.sh <cell-group>                                                        ||
## - execute from pipeline's main  directory                                                                          ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks at cell group level, on all samples that belong to the same cell type using   ||
## MACS3 in paired-end mode.                                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <cell-group> Cell group of samples to call peaks on                                                          ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.sorted.broadPeak.filt                                                      ||
## - output will be peaks called in all samples that belong to the input cell type (that passed stage 1 and 2 of QC.  ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE: <sampleName>.filt.nodup.bam                               ||
## - MACS3 and bedtools installed in a conda environment                                                              ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, BLACKLIST                                                       ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##


mapfile -t SAMPLES < ${META_DIR}/downAnalysisSamples_${GROUP}.txt

MACS_GROUP="${macs_group_thres}"
OUTPUT_DIR="${PEAK_DIR_CC}"
PEAKS_NAME="${GROUP}_${MACS_GROUP}"

bamFiles=()

for sampleName in "${SAMPLES[@]}"; do
  if [[ "$sampleName" =~ _merged$ ]]; then
    bamFiles+=("${ALIGNED_DIR}/${sampleName}.sorted.bam")
  else
    bamFiles+=("${ALIGNED_DIR}/${sampleName}.filt.nodup.bam")
  fi
done

echo ${#SAMPLES[@]}

INPUTDIR="/lustre/projects/Research_Project-MRC190311/ATACSeq/methodsPaper/3_aligned/merged/passQC/"
mergedBam="${INPUTDIR}/${GROUP}.sorted.merged.bam"

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

echo "Starting peak calling using MACS3 PE on samples that belong to cell fraction ${GROUP} at: "
date -u
echo "Merged bam file is: ${mergedBam}"

echo "Q-value threshold is ${MACS_GROUP}"

cd ${TMPDIR}

#macs3 callpeak -t ${bamFiles[@]} --outdir ${OUTPUT_DIR} -n ${PEAKS_NAME} -f BAMPE -g 2.9e9 -q ${MACS_GROUP} --keep-dup all --nomodel

macs3 callpeak -t ${mergedBam} --outdir ${OUTPUT_DIR} -n ${PEAKS_NAME} -f BAMPE -g 2.9e9 -q ${MACS_GROUP} --keep-dup all --nomodel

echo "Number of peaks before filtering: " 
wc -l "${OUTPUT_DIR}/${PEAKS_NAME}_peaks.narrowPeak"

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a "${OUTPUT_DIR}/${PEAKS_NAME}_peaks.narrowPeak" -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | awk '!/^(chrY|chrX)/' > "${OUTPUT_DIR}/${PEAKS_NAME}.narrowPeak.filt"

echo "Number of peaks after filtering: " 
wc -l "${OUTPUT_DIR}/${PEAKS_NAME}.narrowPeak.filt"

echo "Sorting peaks by chr"
sort -k1 "${OUTPUT_DIR}/${PEAKS_NAME}.narrowPeak.filt" > "${OUTPUT_DIR}/${PEAKS_NAME}.sorted.narrowPeak.filt"

if [[ ! -f "${OUTPUT_DIR}/${PEAKS_NAME}.narrowPeak.filt" ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 9.0 was properly run." ; exit 1; }
else
  echo "Finished calling peaks using MACS3 PE on samples in group ${GROUP}"
  echo Job ended on:
  date -u
fi