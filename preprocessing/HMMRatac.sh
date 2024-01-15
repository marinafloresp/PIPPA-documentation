## Performs sample level peak calling using HMMRATAC with paired end BAM files


## EXECUTION
# sh ./ATACSeq/preprocessing/HMMRatac.sh <sampleName>
# where 
# <sampleName> is the sample specific file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST, HMMRATAC

## REQUIRES the following software
# samtools

## INPUT
# ${sampleName}.sorted.bam.bai
# ${sampleName}.filt.nodup.bam

## OUTPUT
# ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak


cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}/HMMRATAC
sampleName=$1

alignSample=${sampleName}.filt.nodup.bam
baiSample=${sampleName}.filt.nodup.bam.bai
blackListBED="/lustre/projects/Research_Project-MRC190311/references/ATAQSeq/hg38.blacklist.bed"

software="/lustre/projects/Research_Project-MRC190311/software"

cd $software
echo "Current directory:" $software
#samtools view -H ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' >  #${PEAKDIR}/HMMRATAC/genome.info

java -Xmx4G -jar HMMRATAC_V1.2.10_exe.jar -b ${ALIGNEDDIR}/${alignSample} -i ${ALIGNEDDIR}/${baiSample} -g ${GENOMESIZE} -o ${PEAKDIR}/HMMRATAC/${sampleName} -e $blackListBED

awk -v OFS="\t" '$13>=10 {print}' ${PEAKDIR}/HMMRATAC/${sampleName}_peaks.gappedPeak > ${PEAKDIR}/HMMRATAC/${sampleName}.filteredPeaks.gappedPeak

echo "Peaks called"