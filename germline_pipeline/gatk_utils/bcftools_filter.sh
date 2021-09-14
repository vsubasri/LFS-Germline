#!/bin/bash
#PBS -S /bin/bash
#PBS -N bcftools-Filtering
#PBS -o log/bcftools-filtering.out
#PBS -e log/bcftools-filtering.err
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 bcftools/1.6

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar
TMPDIR=/localhd/$PBS_JOBID

echo $sample_path

echo "Filtering sample $sample using bcftools"

bcftools view $sample_path \
        | bcftools norm -m - \
        | bcftools norm -f $ref_fasta \
        | bcftools filter -e 'AD[1-] < 5' \
        | bcftools filter -i 'TYPE!="snp"' \
	| bcftools filter -e 'FORMAT/GT=="./."' \
	| bcftools filter -i 'AF[1] > 0.05' \
        | bcftools view -o ${indel_dir}${sample}_indel.vcf

bcftools view $sample_path \
        | bcftools norm -m - \
        | bcftools norm -f $ref_fasta \
        | bcftools filter -e 'AD[1-] < 5' \
        | bcftools filter -i 'TYPE="snp"' \
	| bcftools filter -e 'FORMAT/GT=="./."' \
	| bcftools view -o ${snv_dir}${sample}_SNV.vcf
