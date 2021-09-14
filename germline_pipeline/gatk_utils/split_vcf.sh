#!/bin/bash
#PBS -S /bin/bash
#PBS -N splitvcf
#PBS -e log/splitvcf.err
#PBS -o log/splitvcf.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bcftools/1.6

IFS= read -a array <<< $(grep "#CHROM" $multisample | head -1 | awk '{for(i=10;i<=NF;++i)print $i}')
samples=${array[0]}
for i in `seq 1 $(echo "$samples" | wc -w)`;
do
	sample=$(echo "$samples" | cut -d" " -f$i)
	sample_path=${gatkdir}/${sample}/${sample}.vcf
	bcftools view -s $sample $multisample > $sample_path
#	bcftools view -s $(echo "$samples" | cut -d" " -f$i) $multisample > $gatkdir"/"$(echo "$samples" | cut -d" " -f$i)"/"$(echo "$samples" | cut -d" " -f$i).vcf
done
