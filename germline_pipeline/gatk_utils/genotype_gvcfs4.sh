#!/bin/bash
#PBS -S /bin/bash
#PBS -N GenotypeGVCFs
#PBS -e log/genotypeGVCFs.err
#PBS -o log/genotypeGVCFs.out
#PBS -l mem=100g
#PBS -l vmem=100g
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91

####################### GATK CombineGVCFs   ##################################################################################
###
###
# Take as input the per sample per chromosome gVCF files and produces the combined gvcf file which is
# meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs
# Reference: https://gatkforums.broadinstitute.org/gatk/discussion/10061/using-genomicsdbimport-to-consolidate-gvcfs-for-input-to-genotypegvcfs-in-gatk4
###
###

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar
TMPDIR=/localhd/$PBS_JOBID

BATCH=$(echo $(basename $INPUTLIST) | sed 's/input.batch\(.*\).list/\1/') 
GATKDIR=$(dirname $INPUTLIST)/

cd $GATKDIR

echo "Batch: $BATCH"
echo "Working Directory: $GATKDIR"

function CombineGVCFs {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx84G $gatk CombineGVCFs \
        --variant ${GATKDIR}input.batch${BATCH}.list \
        -R ${ref_fasta} \
        -O ${GATKDIR}combined_variants.${BATCH}.g.vcf
}
function CombineGVCFs_AS {
     java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk CombineGVCFs \
        -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \ 
        --variant ${GATKDIR}input.batch${BATCH}.list \
        -R ${ref_fasta} \
        -O ${GATKDIR}combined_variants.${BATCH}.g.vcf
}

function GenotypeGVCFs {
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${GATKDIR}combined_variants.${BATCH}.g.vcf \
      -O ${GATKDIR}genotyped_variants.${BATCH}.g.vcf \
      --dbsnp ${db_snp} \
      --annotation InbreedingCoeff \
      --annotation QualByDepth \
      --annotation MappingQualityRankSumTest \
      --annotation ReadPosRankSumTest \
      --annotation FisherStrand
}

function GenotypeGVCFs_AS {
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${GATKDIR}combined_variants.${BATCH}.g.vcf \
      -O ${GATKDIR}genotyped_variants.${BATCH}.g.vcf \
      --dbsnp ${db_snp} \
      -G StandardAnnotation \
      -G AS_StandardAnnotation \
      -G StandardHCAnnotation \
      --annotation InbreedingCoeff \
      --annotation QualByDepth \
      --annotation MappingQualityRankSumTest \
      --annotation ReadPosRankSumTest \
      --annotation FisherStrand
}


function GenotypeGVCFs_intlist {
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${GATKDIR}combined_variants.${BATCH}.g.vcf \
      -O ${GATKDIR}genotyped_variants.${BATCH}.g.vcf \
      -L ${interval_list} \
      --annotation InbreedingCoeff \
      --annotation QualByDepth \
      --annotation MappingQualityRankSumTest \
      --annotation ReadPosRankSumTest \
      --annotation FisherStrand \
      --dbsnp ${db_snp} \
      -ip 100
}


function GenotypeGVCFs_intlist_AS {
  java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${GATKDIR}combined_variants.${BATCH}.g.vcf \
      -O ${GATKDIR}genotyped_variants.${BATCH}.g.vcf \
      -L ${interval_list} \
      --dbsnp ${db_snp} \
      -G StandardAnnotation \
      -G AS_StandardAnnotation \
      -G StandardHCAnnotation \
      --annotation InbreedingCoeff \
      --annotation QualByDepth \
      --annotation MappingQualityRankSumTest \
      --annotation ReadPosRankSumTest \
      --annotation FisherStrand \
      -ip 100    
}

################################################## SANITY CHECK ##################################################

echo "
        
	il=${INTERVALLIST}
        as=${ALLELESPECIFIC}
        gatk=${GATKDIR}
"

##################################################################################################################

if ([ "${INTERVALLIST}" == "True" ] &&  [ "${ALLELESPECIFIC}" == "True" ]); then
	CombineGVCFs_AS
	GenotypeGVCFs_intlist_AS
elif ([ "${INTERVALLIST}" == "True" ] &&  [ "${ALLELESPECIFIC}" != "True" ]); then
	CombineGVCFs
	GenotypeGVCFs_intlist
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" == "True" ]); then
	CombineGVCFs_AS
	GenotypeGVCFs_AS
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" != "True" ]); then
	CombineGVCFs
	GenotypeGVCFs
else
   echo "Error in Genotyping"
   exit
fi

exit_status=$?;
echo EXIT STATUS: ${exit_status}
echo END: `date`
exit ${exit_status};
