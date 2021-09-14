#!/bin/bash
#PBS -S /bin/bash
#PBS -N AnnotateCNVs
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load python/3.5.2

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/processCNV.py -a /hpf/largeprojects/adam/projects/lfs/israeli/cnvs/
