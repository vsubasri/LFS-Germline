#!/bin/bash

input_dir=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/kics.run.list
gatk_dir=/hpf/largeprojects/adam/projects/kics/kics_germline/gatk

if [[ -d $input_dir ]]; then
	files=( $(find $input_dir -type f -name "*.bam") )
else
	mapfile -t files < $input_dir
fi

for bam in ${files[@]}
do
	sample=$(basename $bam .realigned-recalibrated.bam)
	sample_vcf=$gatk_dir/$sample/$sample.vcf
	if [ ! -f $sample_vcf ] ; then
		echo $sample_vcf
	fi
done
 
