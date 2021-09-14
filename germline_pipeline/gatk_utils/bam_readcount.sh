#!/bin/bash
#PBS -S /bin/bash
#PBS -N bam-readcount
#PBS -e bam-readcount.err
#PBS -o bam-readcount.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bam-readcount/0.8.0

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/gatk/filtered_variants/
site_list=
#awk -F '\t' -v OFS='\t' '{print $1, $2, $2}' filtered_snvs.tab > snvs_site_list
#tail -n +2 snvs_site_list

bam-readcount -f $ref_fasta -q 10 -b 15 -w 1 -l $site_list $bam > ${outdir}${sample}.readcount
