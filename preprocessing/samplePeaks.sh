## Performs sample level peak calling using MACS with 1. shifted tagAlign files 
## and 2. paired end BAM files
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ATACSeq/preprocessing/5_samplePeaks.sh <sampleName>
# where 
# <sampleName> is the sample specific file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# shifted tag align file
# .filt.nodup.bam

## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/*.broadPeak (and other macs output)
# ${PEAKDIR}/MACS/BAMPE/*.broadPeak
# ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt
# ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt

cd ${ALIGNEDDIR}
sampleName=$1

## ==========================
## MACS2 TA mode peak calling
## ==========================

mkdir -p ${PEAKDIR}/MACS/ShiftedTagAlign
f_TA=${sampleName}.tn5.tagAlign.gz
#if [ ! -s ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt ]	
##then
echo "Calling peaks using MACS2 TA for"" ${sampleName}"
date -u	
macs2 callpeak -t ${f_TA} --outdir ${PEAKDIR}/MACS/ShiftedTagAlign -n ${sampleName} -f BED -g 2.9e9 -q 5e-2 --keep-dup all --shift 100 --extsize 200 --nomodel --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
| grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt

rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak
rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.xls
rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.gappedPeak  

#else
#	{ echo "TA Filtered peak file found so not calling peaks"; exit 1;}
#fi

## ============================= 
## MACS2 BAMPE mode peak calling
## ============================= 
mkdir -p ${PEAKDIR}/MACS/BAMPE
f_PE=${sampleName}.filt.nodup.bam
#if [ ! -s ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]	
#then
echo "Calling peaks using MACS2 BAMPE for"" ${sampleName}"
date -u	
macs2 callpeak -t ${f_PE} --outdir ${PEAKDIR}/MACS/BAMPE -n ${sampleName} -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
| grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt

rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.xls
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.gappedPeak  
  
#else
#	{ echo "BAMPE Filtered peak file found so not calling peaks"; exit 1;}
#fi
