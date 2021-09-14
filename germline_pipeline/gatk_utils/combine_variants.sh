#!/bin/bash
#PBS -S /bin/bash
#PBS -N CombineVariants
#PBS -e log/combinevariants.err
#PBS -o log/combinevariants.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 picard-tools/2.5.0

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
gatk3=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=/localhd/$PBS_JOBID
final_dir=${var_dir}final_variants/
fpfilter_dir=${var_dir}fpfilter_snvs/
indel_dir=${var_dir}sample_indels/
outdir=${var_dir}multisample/

mkdir $outdir

find ${fpfilter_dir} -name "*fpfilter.vcf" > ${fpfilter_dir}input.list

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk3 \
   -T CombineVariants \
   -R $ref_fasta \
   --variant ${fpfilter_dir}input.list \
   -o ${outdir}fpfiltered_all.vcf \
   -genotypeMergeOptions UNIQUIFY

find ${indel_dir} -name "*indel.vcf" > ${indel_dir}input.list

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk3 \
   -T CombineVariants \
   -R $ref_fasta \
   --variant ${indel_dir}input.list \
   -o ${outdir}indel_all.vcf \
   -genotypeMergeOptions UNIQUIFY

find ${final_dir} -name "*.vcf" > ${final_dir}input.list

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk3 \
   -T CombineVariants \
   -R $ref_fasta \
   --variant ${final_dir}input.list \
   -o ${outdir}final_variants.vcf \
   -genotypeMergeOptions UNIQUIFY

