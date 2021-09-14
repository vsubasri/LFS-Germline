#!/bin/bash

#PBS -S /bin/bash
#PBS -N HaplotypeCaller
#PBS -e log/haplotypecaller.err
#PBS -o log/haplotypecaller.out
#PBS -l mem=32g
#PBS -l vmem=48g
#PBS -l walltime=300:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar
TMPDIR=/localhd/$PBS_JOBID

sample=$(basename $PBS_O_INITDIR)
sample_dir=$PBS_O_INITDIR
echo $sample_dir
cd $sample_dir

function hc {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk HaplotypeCaller \
	-R $ref_fasta \
	-I $bam \
	-O $sample_dir/$sample.g.vcf \
	--base-quality-score-threshold 20 \
	--emit-ref-confidence GVCF \
	-D ${db_snp} 
}

function hc_intlist {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk HaplotypeCaller \
        -R $ref_fasta \
        -I $bam \
        -O $sample_dir/$sample.g.vcf \
        --base-quality-score-threshold 20 \
        --emit-ref-confidence GVCF \
        -L ${interval_list} \
        -ip 100 \
        -D ${db_snp}
}

function hc_intlist_allelespecific {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk HaplotypeCaller \
        -R $ref_fasta \
        -I $bam \
        -O $sample_dir/$sample.g.vcf \
        --base-quality-score-threshold 20 \
        --emit-ref-confidence GVCF \
        -L ${interval_list} \
        -ip 100 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -D ${db_snp}
}

function hc_allelespecific {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk HaplotypeCaller \
        -R $ref_fasta \
        -I $bam \
        -O $sample_dir/$sample.g.vcf \
        --base-quality-score-threshold 20 \
        --emit-ref-confidence GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -D ${db_snp}
}

################################################## SANITY CHECK ##################################################

echo "

        il=${INTERVALLIST}
        as=${ALLELESPECIFIC}
        bam=${bam}
"

##################################################################################################################

if ([ "${INTERVALLIST}" == "True" ] && [ "${ALLELESPECIFIC}" == "True" ]); then
	hc_intlist_AS
elif ([ "${INTERVALLIST}" == "True" ] && [ "${ALLELESPECIFIC}" != "True" ]); then
	hc_intlist
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" == "True" ]); then
	hc_AS
elif ([ "${INTERVALLIST}" != "True" ] &&  [ "${ALLELESPECIFIC}" != "True" ]); then
	hc
else 
  echo "Error in Genotyping"
  exit
fi

exit_status=$?;
echo EXIT STATUS: ${exit_status}
echo END: `date`
exit ${exit_status};
