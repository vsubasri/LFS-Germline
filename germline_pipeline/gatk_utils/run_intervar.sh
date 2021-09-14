#!/bin/bash
#PBS -S /bin/bash
#PBS -N InterVar
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -l walltime=248:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load annovar/2018.04.17 python/3.5.2

cd /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/InterVar/

if [ -f ${output}.hg19_multianno.txt ] && [ ! -f ${output}.hg19_multianno.txt.intervar ]; then
	echo "${output}.hg19_multianno.txt already exists, skipping ANNOVAR"
	/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/InterVar/Intervar.py --skip_annovar -b hg19 -i ${input} --input_type=VCF  -o ${output} -d /hpf/largeprojects/adam/projects/lfs/resources/humandb
elif [ -f ${output}.hg19_multianno.txt.intervar ] ; then
	echo "${output}.hg19_multianno.txt.intervar already exists, skipping sample"
else
	echo "Annotating sample: $output"
        /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/InterVar/Intervar.py -b hg19 -i ${input} --input_type=VCF  -o ${output} -d /hpf/largeprojects/adam/projects/lfs/resources/humandb
fi
