#!/bin/bash
#PBS -S /bin/bash
#PBS -N snpeff
#PBS -e snpeff.err
#PBS -o snpeff.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load java/1.8.0_91

TMPDIR=/localhd/$PBS_JOBID
config=/hpf/largeprojects/adam/projects/lfs/resources/snpEff/snpEff.config
snpeff=/hpf/largeprojects/adam/projects/lfs/resources/snpEff/snpEff.jar
input=/hpf/largeprojects/adam/projects/lfs/lfs_germline/gatk/filtered_variants/het_hom_genotypes/fpfilter_snvs/fpfiltered_snvs.vcf
ens_output=/hpf/largeprojects/adam/projects/lfs/lfs_germline/gatk/filtered_variants/het_hom_genotypes/fpfilter_snvs/fpfiltered_snvs.refseq_ens_annotated.vcf
refseq_output=/hpf/largeprojects/adam/projects/lfs/lfs_germline/gatk/filtered_variants/het_hom_genotypes/fpfilter_snvs/fpfiltered_snvs.refseq_annotated.vcf

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G ${snpeff} eff \
	-v -i vcf -o vcf -c ${config} \
	-spliceSiteSize 7 hg19 \
	${input} > ${refseq_output}

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G ${snpeff} eff \
	-v -motif -nextprot -i vcf -o vcf -c ${config} GRCh37.87 \
	${refseq_output} > ${ens_output}	



