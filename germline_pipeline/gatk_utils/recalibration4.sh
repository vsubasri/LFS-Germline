#!/bin/bash
#PBS -S /bin/bash
#PBS -N Recalibration
#PBS -e log/VQSR.err
#PBS -o log/VQSR.out
#PBS -l mem=32g
#PBS -l vmem=100g
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
hapmap=/hpf/largeprojects/adam/projects/lfs/resources/hapmap_3.3.b37.vcf
omni=/hpf/largeprojects/adam/projects/lfs/resources/1000G/1000G_omni2.5.b37.vcf
thouG=/hpf/largeprojects/adam/projects/lfs/resources/1000G/1000G_phase1.snps.high_confidence.b37/1000G_phase1.snps.high_confidence.b37.vcf
mills=/hpf/largeprojects/adam/projects/lfs/resources/Mills_and_1000G_gold_standard.indels.b37.vcf
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar
TMPDIR=/localhd/$PBS_JOBID

BATCH=$(echo $(basename $INPUTFILE) | sed 's/genotyped_variants.\(.*\).g.vcf/\1/')
GATKDIR=$(dirname $INPUTFILE)/

cd $GATKDIR

####################### VQSR ##################################################################################

# Methods Reference: https://software.broadinstitute.org/gatk/documentation/article.php?id=39
# Parameter Reference: https://software.broadinstitute.org/gatk/documentation/article.php?id=1259
# Tutorial: https://software.broadinstitute.org/gatk/documentation/article.php?id=2805

function snp_VQSR {
  ##Step 1: VariantRecalibrator
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -R ${ref_fasta} \
    -V ${INPUTFILE} \
    -O ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
    --resource omni,known=false,training=true,truth=true,prior=12.0:${omni} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:${thouG} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
    -mode SNP \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -R ${ref_fasta} \
    -O ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf  \
    -V ${INPUTFILE} \
    --recal-file ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP
}

function indel_VQSR {
  ##Step 1: VariantRecalibrator

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -R ${ref_fasta} \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
    -mode INDEL \
    -O ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -R ${ref_fasta} \
    -O ${GATKDIR}indel.recalibrated.${BATCH}.vcf \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    --recal-file ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL
}

function snp_VQSR_AS {
  ##Step 1: VariantRecalibrator
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -AS \
    -an AS_QD -an AS_FS -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an AS_InbreedingCoeff \
    -R ${ref_fasta} \
    -V ${INPUTFILE} \
    -O ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
    --resource omni,known=false,training=true,truth=true,prior=12.0:${omni} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:${thouG} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -mode SNP \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.SNPs.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -AS \
    -R ${ref_fasta} \
    -O ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf  \
    -V ${INPUTFILE} \
    --recal-file ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP
}

function indel_VQSR_AS {
  ##Step 1: VariantRecalibrator

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -AS \
    -an AS_QD -an AS_FS -an AS_ReadPosRankSum -an AS_MQRankSum -an AS_SOR -an AS_InbreedingCoeff \
    -R ${ref_fasta} \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -O ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -mode INDEL \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.indels.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -AS \
    -R ${ref_fasta} \
    -O ${GATKDIR}indel.recalibrated.${BATCH}.vcf \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    --recal-file ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL
}

function snp_VQSR_intlist {
  ##Step 1: VariantRecalibrator
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -R ${ref_fasta} \
    -V ${INPUTFILE} \
    -O ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    -L ${interval_list} \
    -ip 100 \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
    --resource omni,known=false,training=true,truth=true,prior=12.0:${omni} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:${thouG} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
    -mode SNP \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -R ${ref_fasta} \
    -O ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -V ${INPUTFILE} \
    -L ${interval_list} \
    -ip 100 \
    --recal-file ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP
}

function indel_VQSR_intlist {
  ##Step 1: VariantRecalibrator

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -R ${ref_fasta} \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -O ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    -L ${interval_list} \
    -ip 100 \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
    -mode INDEL \
    -tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    -rscript-file ${GATKDIR}recal.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -R ${ref_fasta} \
    -O ${GATKDIR}indel.recalibrated.${BATCH}.vcf \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -L ${interval_list} \
    -ip 100 \
    --recal-file ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL
}

function snp_VQSR_AS_intlist {
  ##Step 1: VariantRecalibrator
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -AS \
    -an AS_QD -an AS_FS -an AS_ReadPosRankSum -an AS_MQ -an AS_MQRankSum -an AS_SOR -an AS_InbreedingCoeff \
    -R ${ref_fasta} \
    -V ${INPUTFILE} \
    -L ${interval_list} \
    -ip 100 \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
    --resource omni,known=false,training=true,truth=true,prior=12.0:${omni} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:${thouG} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -mode SNP \
    -O ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.SNPs.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -AS \
    -R ${ref_fasta} \
    -O ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf  \
    -V ${INPUTFILE} \
    -L ${interval_list} \
    -ip 100 \
    --recal-file ${GATKDIR}raw.SNPs.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.SNPs.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP
}

function indel_VQSR_AS_intlist {
  ##Step 1: VariantRecalibrator

  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantRecalibrator \
    -AS \
    -an AS_QD -an AS_FS -an AS_ReadPosRankSum -an AS_MQRankSum -an AS_SOR -an AS_InbreedingCoeff \
    -R ${ref_fasta} \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -L ${interval_list} \
    -ip 100 \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${db_snp} \
    -mode INDEL \
    -O ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --rscript-file ${GATKDIR}recal.indels.plots.${BATCH}.R

  ##Step 2: ApplyRecalibration
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk ApplyVQSR \
    -AS \
    -R ${ref_fasta} \
    -O ${GATKDIR}indel.recalibrated.${BATCH}.vcf \
    -V ${GATKDIR}SNPs.recalibrated.${BATCH}.vcf \
    -L ${interval_list} \
    -ip 100 \
    --recal-file ${GATKDIR}raw.indels.recal.${BATCH}.vcf \
    --tranches-file ${GATKDIR}raw.indels.${BATCH}.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL
}

################################################## SANITY CHECK ##################################################

echo "

        il=${INTERVALLIST}
        as=${ALLELESPECIFIC}
        gatk=${GATKDIR}

"

##################################################################################################################

if ([ "${INTERVALLIST}" == "True" ] &&  [ "${ALLELESPECIFIC}" == "True" ]); then
	snp_VQSR_AS_intlist
	indel_VQSR_AS_intlist
elif ([ "${INTERVALLIST}" == "True" ] &&  [ "${ALLELESPECIFIC}" != "True" ]); then
	snp_VQSR_intlist
	indel_VQSR_intlist
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" == "True" ]); then
	snp_VQSR_AS
	indel_VQSR_AS
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" != "True" ]); then
	snp_VQSR
	indel_VQSR
else 
	echo "Error in Recalibration"
	exit
fi

exit_status=$?;
echo EXIT STATUS: ${exit_status}
echo END: `date`
exit ${exit_status};
