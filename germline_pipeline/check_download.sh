#!/bin/bash

input_dir=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/all.files
gatk_dir=/hpf/largeprojects/adam/projects/kics/kics_germline/gatk

if [[ -d $input_dir ]]; then
	files=( $(find $input_dir -type f -name "*.bam.md5") )
else
	mapfile -t files < $input_dir
fi

for bammd in ${files[@]}
do
	sample=$(basename $bammd .realigned-recalibrated.bam.md5)
	sample_dir=$gatk_dir/$sample
	if [ ! -d $sample_dir ] ; then
		bam=$(echo $bammd | cut -f1 -d".")
		tombstone download ${bam}.realigned-recalibrated.bam.tombstone
	fi
done
 
