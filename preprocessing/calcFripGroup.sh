## calculates number of reads within peaks for group level peaks

## EXECUTION
# sh ./ATACSeq/preprocessing/calcFripGroup.sh 
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools, samtools

## INPUT
# samples in group

## OUTPUT
# ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv

FAIL=$@

samples=("${PEAKDIR}/MACS/ShiftedTagAlign/subset/*.broadPeak.filt")
echo $samples

for sampleName in ${FAIL[@]}
do
  ##Formulate file name
  tagAlign=${sampleName}.tn5.tagAlign.gz

	## save output in txt file
	echo "sampleName,tagAlignTotalReads,sampleTagAlignPeaks,readsInSampleTagAlignPeaks,subsetTagAlignPeaks,readsInSubsetTagAlignPeaks" > ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#sampletotalreads
	echo -n $(zcat ${ALIGNEDDIR}/${tagAlign} | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#samplepeaks
	echo -n $(wc -l ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
  #reads in peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${tagAlign} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#subsetpeaks
	echo -n $(wc -l ${PEAKDIR}/MACS/ShiftedTagAlign/subset/*.broadPeak | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#reads in subset peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/MACS/ShiftedTagAlign/subset/*.broadPeak | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${tagAlign} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
done

