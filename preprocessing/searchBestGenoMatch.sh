


sampleName=$1
sampleName=${sampleName%.selfSM}
vcfid=$2
vcfid="${vcfid%"${vcfid##*[![:space:]]}"}" 

cd ${ALIGNEDDIR}/

echo "processing" ${sampleName} "with" ${vcfid}
#"${GENODIR}/hg38/allchr_filt_rsq_maf.vcf.gz"
${VERIFYBAMID} --vcf "/lustre/projects/Research_Project-MRC190311/SNPArray/bdrBrain/Marina16SamplesOnly_vcf.vcf" --bam baseRecalibrate/${sampleName}_baserecal.bam --out genotypeConcordance/${sampleName} --ignoreRG --smID ${vcfid} --best

