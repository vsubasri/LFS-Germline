#!/bin/bash
#PBS -S /bin/bash
#PBS -N PostFiltering
#PBS -l mem=10g
#PBS -l vmem=10g
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load python/3.5.2 bedtools/2.27.1

python /hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/structural_variations/denovo_variants.py


