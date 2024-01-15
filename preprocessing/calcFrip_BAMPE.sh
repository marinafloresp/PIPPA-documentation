## calculates number of reads within peaks for sample level peaks

## EXECUTION
# sh ./ATACSeq/preprocessing/6_calcFrip.sh <sampleName>
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
# ${PEAKDIR}/QCOutput/FRIP_*.csv

sampleName=$1

## save output in txt file
echo "SampleName,BAMTotalReads,MACS2PEPeaks,ReadsinMACS2PEPeaks,FRIP" > ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv

echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv

# formulate peak input filename
bam=${sampleName}.filt.nodup.bam 

if [ -s ${ALIGNEDDIR}/${bam} ]
then 
  total_reads=$(samtools view -c ${ALIGNEDDIR}/${bam})
	echo -n $total_reads, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
else
	echo -n NA, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
fi

# count reads in peak regions

## MACS2 peaks called from bam files with paired end flag

if [ -s ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]
then
	echo -n $(wc -l ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
  peaks=$(bedtools sort -i ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${bam} -b stdin -ubam | samtools view -c)
	echo -n $peaks, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
  fraction=$(echo "scale=4; $peaks / $total_reads" | bc -l)
  echo -n "$fraction", >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
else
	echo NA,NA,NA, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
fi
