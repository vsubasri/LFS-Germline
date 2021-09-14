#!/bin/bash
#PBS -S /bin/bash
#PBS -N PostFiltering
#PBS -e log/PostFiltering.err
#PBS -o log/PostFiltering.out
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load python/3.5.2

#intervar_dir=/hpf/largeprojects/adam/projects/kics/kics_germline/gatk/filtered_variants/intervar
intervar_dir=/hpf/largeprojects/adam/projects/lfs/resources/1000G/1000G_phase3_v5_20130502/chr15

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/post-filtering.py -i $intervar_dir 

