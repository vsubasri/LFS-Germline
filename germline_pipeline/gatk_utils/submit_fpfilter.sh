#!/bin/bash
#PBS -S /bin/bash
#PBS -N fpfilter-Filtering
#PBS -o log/fpfilter.out
#PBS -e log/fpfilter.err
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l walltime=64:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load perl/5.20.1 bam-readcount/0.8.0

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
fpfilter=/hpf/largeprojects/adam/projects/lfs/code/gatk/post-analysis/filter-annot-pipeline/fpfilter.v2.pl

#perl ${fpfilter} --var-file ${vcf} --readcount-file ${brc} --output-file ${output} 
echo "Filtering sample $sample using fpfilter..."
echo "BAM: $bam"
echo "BCFTOOLS: $vcf"
echo "FPFILTER: $output"

perl ${fpfilter} --vcf-file ${vcf} --bam-file ${bam} --sample ${sample} --reference ${ref_fasta} --output ${output}
