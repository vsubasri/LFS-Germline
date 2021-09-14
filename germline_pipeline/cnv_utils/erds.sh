#!/bin/bash
#PBS -S /bin/bash
#PBS -N ERDS
#PBS -e log.erds.err
#PBS -o log.erds.out
#PBS -l vmem=48g
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load erds/1.1

genome=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
erds=/hpf/largeprojects/adam/projects/lfs/resources/erds1.1/erds_pipeline.pl

perl ${erds} -o ${erds_output} -b ${bam_file} -v ${vcf} -r ${genome}
