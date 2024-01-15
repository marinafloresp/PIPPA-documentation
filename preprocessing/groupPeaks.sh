## Performs fraction/group level peak calling using MACS with shifted tagAlign files 
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ATACSeq/preprocessing/groupPeaks.sh 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST, PEAKGROUP

## REQUIRES the following software
# bedtools, macs2

## INPUT
# array of shifted tag align files

## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/subset/${GROUP}_atac.broadPeak.filt (and other macs output)

cd ${ALIGNEDDIR}

FILES=$@
#For Shifted TA mode
#FILES=( "${QC1SAMPLES[@]/%/.tn5.tagAlign.gz}" )

#For BAMPE mode
#FILES=("${QC1SAMPLES[@]/%/.filt.nodup.bam}")

#SAMPLEFILE=$1
#QC1SAMPLES=($(awk '{print $1}' ${SAMPLEFILE} ))
#SUBSET=("${QC1SAMPLES[@]/%/.tn5.tagAlign.gz}" )


#Peaks called using Tag Align samples
#echo 'Files for peak-calling are:' ${FILES[@]}
#macs2 callpeak -t ${FILES[@]} --outdir ${PEAKGROUP} -n ${GROUP}_TA_atac -f BED -g 2.9e9 -q 5e-2 --keep-dup all  --shift 100 --extsize 147 --nomodel --broad --broad-cutoff 5e-2

#Peaks called using PE samples
echo 'Files for peak-calling are:' ${FILES[@]}
macs2 callpeak -t ${FILES[@]} --outdir ${PEAKDIR}/MACS/BAMPE/subset -n ${GROUP}_BAMPE_atac -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --nomodel --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/subset/${GROUP}_BAMPE_atac_peaks.broadPeak -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	| grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/BAMPE/subset/${GROUP}_BAMPE_atac.broadPeak.filt

echo 'Finished calling peaks in common'
