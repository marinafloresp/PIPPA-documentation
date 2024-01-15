## calculates number of reads within peaks for sample level peaks using HMMRATAC called peaks

## EXECUTION
# sh ./ATACSeq/preprocessing/calcFrip_HMMR.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools, samtools

## INPUT
# .filt.nodup.bam

## OUTPUT
# ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_*.csv


sampleName=$1

## save output in txt file
echo "SampleName,BAMTotalReads,HMMRATACPeaks,ReadsinHMMRATACPeaks,FRIP" > ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv

echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv

# formulate peak input filename
bam=${sampleName}.filt.nodup.bam 

if [ -s ${ALIGNEDDIR}/${bam} ]
then 
  total_reads=$(samtools view -c ${ALIGNEDDIR}/${bam})
	echo -n $total_reads, >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
else
	echo -n NA, >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
fi

## HMMRATAC peaks called from bam files 

if [ -s ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak ]
then
	echo -n $(wc -l ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
  peaks=$(bedtools sort -i ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${bam} -b stdin -ubam | samtools view -c)
	echo -n $peaks, >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
  fraction=$(echo "scale=4; $peaks / $total_reads" | bc -l)
  echo -n "$fraction", >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
else
	echo NA,NA,NA, >> ${PEAKDIR}/QCOutput/HMMRATAC/FRIP_${sampleName}.csv
fi
