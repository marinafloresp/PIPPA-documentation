## calculates number of reads within peaks for group level peaks called by HMMRATAC

## EXECUTION
# sh ./ATACSeq/preprocessing/calcFripGroup_HMMR.sh 
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
# ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}.csv

SAMPLES=$@
echo $sample

#HMMRATAC peaks files
peaks=("${PEAKDIR}/HMMRATAC/subset/*.filteredPeaks.gappedPeak")


for sampleName in ${SAMPLES[@]}
do
	sample=${sampleName}.filt.nodup.bam
  ## save output in txt file
  
  echo "sampleName,HMMRTotalReads,sampleHMMRPeaks,readsInSampleHMMRPeaks, N_samplesPeaks,IRF_samplesPeaks, Sox_samplesPeaks, TN_samplesPeaks, readsInN_samplesPeaks,readsInIRF_samplesPeaks, readsInSox_samplesPeaks, readsInTN_samplesPeaks" > ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}.csv

	echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}.csv
	#sampletotalreads
  total_reads=$(samtools view -c ${ALIGNEDDIR}/${sample})
	echo -n $total_reads, >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}.csv
	#samplepeaks
	echo -n $(wc -l ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}.csv
  #reads in peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak| bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${sample} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}..csv
	#subsetpeaks
  for peak in ${peaks[@]}
  do 
    echo -n $(wc -l $peak | awk '{print $1}'), >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}..csv
  done
	#reads in subset peaks
  for peak in ${peaks[@]}
  do
    frac=$(bedtools sort -i $peak | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${sample} -b stdin | wc -l | cut -f1 -d' ')
    fraction=$(echo "scale=4; $frac / $total_reads" | bc -l)
	  echo -n "$fraction", >> ${PEAKDIR}/QCOutput/subset/FRIP_HMMR_${sampleName}..csv
  done
done