#!/bin/bash
#PBS -S /bin/bash
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -N ApplyFilters
#PBS -e log/ApplyFilters.err
#PBS -o log/ApplyFilters.out

module load python/3.5.2

svdir=/hpf/largeprojects/adam/projects/kics/kics_germline/structural_variations/final_svs

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/applyFilter.py -s ${svdir}

