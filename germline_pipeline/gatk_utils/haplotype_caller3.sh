#!/bin/bash

#PBS -S /bin/bash
#PBS -N HaplotypeCaller
#PBS -e log/haplotypecaller.v3.err
#PBS -o log/haplotypecaller.v3.out
#PBS -V
#PBS -l mem=32g
#PBS -l vmem=100g
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=/localhd/$PBS_JOBID

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T HaplotypeCaller \
	-R $ref_fasta \
	-I $bam \
	-o ${sample_dir}$(basename $bam).g.vcf \
	-L ${interval_list} \
	-ip 100 \
	--emitRefConfidence GVCF \
	--genotyping_mode DISCOVERY \
	-stand_call_conf 30 \
	-rf BadCigar \
	--min_base_quality_score 20 \
	-dfrac 0.99 \
	--dbsnp ${db_snp} 


exit_status=$?;
echo EXIT STATUS: ${exit_status}
echo END: `date`
exit ${exit_status};
