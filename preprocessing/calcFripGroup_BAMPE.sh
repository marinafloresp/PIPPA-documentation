## calculates number of reads within peaks for group level peaks called using BAMPE mode

## EXECUTION
# sh ./ATACSeq/preprocessing/calcFripGroup_BAMPE.sh 
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
# ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv

SAMPLES=$@

#For BAMPE mode 
peaks=("${PEAKDIR}/MACS/BAMPE/subset/*.broadPeak.filt")

for sampleName in ${SAMPLES[@]}
do
  ##Formulate file name
	sample=${sampleName}.filt.nodup.bam
 
  ## save output in txt file
  echo "sampleName,BAMPETotalReads,sampleBAMPEPeaks,readsInSampleBAMPEPeaks, IRF_samplesPeaks,N_samplesPeaks,Sox_samplesPeaks, TN_samplesPeaks, readsInIRF_samplesPeaks, readsInN_samplesPeaks, readsInSox_samplesPeaks, readsInTN_samplesPeaks" > ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv

	echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
	#sampletotalreads
  total_reads=$(samtools view -c ${ALIGNEDDIR}/${sample})
	echo -n $total_reads, >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
	#samplepeaks
	echo -n $(wc -l ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
  #reads in peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${sample} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
	#subsetpeaks
  for peak in ${peaks[@]}
  do 
    echo -n $(wc -l $peak | awk '{print $1}'), >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
  done
	#reads in subset peaks
  for peak in ${peaks[@]}
  do
    frac=$(bedtools sort -i $peak | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${sample} -b stdin | wc -l | cut -f1 -d' ')
    fraction=$(echo "scale=4; $frac / $total_reads" | bc -l)
	  echo -n "$fraction", >> ${PEAKDIR}/QCOutput/subset/FRIP_BAMPE_${sampleName}.csv
  done
done