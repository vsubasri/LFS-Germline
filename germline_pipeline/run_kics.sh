#!/bin/bash
#PBS -S /bin/bash
#PBS -l vmem=24g

bamlist=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/kics.bam.list
outdir=/hpf/largeprojects/adam/projects/lfs/kics_germline/

module load python/3.5.2

cd $outdir

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/germline_pipeline.py -b $bamlist -o $outdir --cnv

