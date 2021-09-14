#!/bin/bash

#input_dir=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/all.files
input_dir=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/kics.bam.list
gatk_dir=/hpf/largeprojects/adam/projects/kics/kics_germline/gatk

if [[ -d $input_dir ]]; then
#	files=( $(find $input_dir -type f -name "*.bam.md5") )
	files=( $(find $input_dir -type f -name "*.bam") )
else
	mapfile -t files < $input_dir
fi

for bammd in ${files[@]}
do
	sample=$(basename $bammd .realigned-recalibrated.bam.md5)
	sample_dir=$gatk_dir/$sample
	if [ ! -d $sample_dir ] ; then
		bam=$(echo $bammd | cut -f1 -d".")
	#	echo ${bam}.realigned-recalibrated.bam
		if [ -f ${bam}.realigned-recalibrated.bam ] ; then
			tombstone free --dry-run ${bam}.realigned-recalibrated.bam.tombstone
		fi
	fi
done
 
