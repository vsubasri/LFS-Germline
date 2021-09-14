#!/bin/bash
#PBS -N CombineSNVsIndels
#PBS -e log/CombineSNVsIndels.err
#PBS -o log/CombineSNVsIndels.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 picard-tools/2.5.0

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
gatk3=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=/localhd/$PBS_JOBID

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk3 \
   -T CombineVariants \
   -R $ref_fasta \
   --variant ${indels} \
   --variant ${snvs} \
   -o ${output} \
   --genotypemergeoption UNSORTED
