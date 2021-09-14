#!/bin/bash

#### INPUT ####
# main = output directory
main=$1
# bam_file = list of bam files or director
bam_file=$2
###############

sample=$(echo $(basename $bam_file) | cut -d '.' -f 1)

if [ ! -d "$main/$sample" ]; then
	mkdir $main/$sample
	cd $main/$sample
#	echo "[ $sample ] Setting up directory"
	cat /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/template/01_configManta.sh|sed -e "s=sample_path=$bam_file=; s/sample_name/$sample/; s=output_dir=$main=" >$main/$sample/01_configManta.sh
	cat /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/template/02_mantasubmit.sh|sed "s/sample_name/$sample/g; s=output_dir=$main=">$main/$sample/02_mantasubmit.sh
fi

if [ ! -f "$main/$sample/results/variants/diploidSV.vcf.gz" ]; then
#	echo "[ $sample ] Running manta"
	dep=$(qsub $main/$sample/01_configManta.sh)
	qsub -W depend=afterok:$dep $main/$sample/02_mantasubmit.sh
fi

