## Performs group level peak calling using HMMRATAC with paired end BAM files


## EXECUTION
# sh ./ATACSeq/preprocessing/groupHMMRatac.sh <sampleName>
# where 
# <sampleName> is the sample specific file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST, HMMRATAC

## REQUIRES the following software
# samtools

## INPUT
# ${sampleName}.filt.nodup.bai
# ${sampleName}.filt.nodup.bam

## OUTPUT
# ${PEAKDIR}/HMMRATAC/


cd ${ALIGNEDDIR}

FILES=$@
#peaks=("${FILES[@]/%/.filt.nodup.bam}")

bamFiles=()
baiFiles=()

for sampleName in ${FILES[@]}
do
  sampleBAM=${sampleName}.filt.nodup.bam
  bamFiles+=( ${sampleBAM})
done
echo ${bamFiles[@]}

#echo 'Files for peak-calling are:' ${FILES[@]}
blackListBED="/lustre/projects/Research_Project-MRC190311/references/ATAQSeq/hg38.blacklist.bed"

software="/lustre/projects/Research_Project-MRC190311/software"

cd ${ALIGNEDDIR}
#samtools merge ${ALIGNEDDIR}/${GROUP}.bam  ${bamFiles[@]}  #-bSo
#samtools sort -@ 4 -o ${ALIGNEDDIR}/${GROUP}.sorted.bam ${ALIGNEDDIR}/${GROUP}.bam
#samtools index ${ALIGNEDDIR}/${GROUP}.sorted.bam 

cd $software
echo "Current directory:" $software
java -jar HMMRATAC_V1.2.10_exe.jar -b ${ALIGNEDDIR}/${GROUP}.sorted.bam -i ${ALIGNEDDIR}/${GROUP}.sorted.bam.bai -g ${GENOMESIZE} -o ${PEAKDIR}/HMMRATAC/${GROUP} -e $blackListBED

awk -v OFS="\t" '$13>=10 {print}' ${PEAKDIR}/HMMRATAC/${GROUP}_peaks.gappedPeak > ${PEAKDIR}/HMMRATAC/${GROUP}.filteredPeaks.gappedPeak

echo "Peaks called"