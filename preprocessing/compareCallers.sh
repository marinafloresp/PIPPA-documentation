## 

## EXECUTION
# sh ./ATACSeq/preprocessing/compareCallers.sh 
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools

## INPUT
# ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv


## OUTPUT
# ${PEAKDIR}/compareCallers_*.csv

SAMPLES=$@
echo $SAMPLES
echo "Difference between MACS2 peaks called from shifted tag align file and MACS2 peaks called from bam files with paired end flag (8th column)"
## save output in txt file
echo "SampleName,TagAlignTotalReads,BAMTotalReads,TAPeaks,ReadsinTAPeaks,BAMPeaks,ReadsinBAMPeaks," > ${PEAKDIR}/QCOutput/compareCallersFRIP.csv

for sampleName in ${SAMPLES[@]}
do
  # formulate peak input filenames
  sampleName=$(basename ${sampleName%.filt.nodup.bam})
  tagAlign=${sampleName}.tn5.tagAlign.gz
  bam=${sampleName}.filt.nodup.bam 
  
  echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  
  # count total reads in peak input files
  if [ -s ${ALIGNEDDIR}/${tagAlign} ]
  then 
  	echo -n $(zcat ${ALIGNEDDIR}/${tagAlign} | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  else 
  	echo -n NA, >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  fi
  if [ -s ${ALIGNEDDIR}/${bam} ]
  then 
  	echo -n $(samtools view -c ${ALIGNEDDIR}/${bam}), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  else
  	echo -n NA, >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  fi  
  
    ## MACS2 peaks called from shifted tag align file
  if [ -s ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt ]
  then
  	echo -n $(wc -l ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  	echo -n $(bedtools sort -i ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${tagAlign} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  else
  	echo -n NA,NA, >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  fi
  
  ## MACS2 peaks called from bam files with paired end flag
  
  if [ -s ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]
  then
  	echo -n $(wc -l ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  	echo $(bedtools sort -i ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${bam} -b stdin -ubam | samtools view -c), >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  else
  	echo NA,NA >> ${PEAKDIR}/QCOutput/compareCallersFRIP.csv
  fi
  
  #echo -n $(awk 'NR >= 2 { print }' < ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv), >> ${PEAKDIR}/QCOutput/compareCallers.csv
  #taPeaks=$(awk -F',' 'NR >= 2 {print $4}' < ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv)
  #baPeaks=$(awk -F',' 'NR >= 2 {print $6}' < ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv)
  #diff="$((taPeaks-baPeaks))"
  #echo "${diff}", >> ${PEAKDIR}/QCOutput/compareCallers.csv
done


