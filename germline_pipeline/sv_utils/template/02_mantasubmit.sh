#!/bin/bash

#PBS -S /bin/bash
#PBS -V
#PBS -l mem=12g
#PBS -l vmem=12g
#PBS -N manta_sample_name
#PBS -o manta_sample_name.log
#PBS -e manta_sample_name.log
#PBS -l nodes=1:ppn=8,walltime=72:00:00

module purge; module load python/2.7.10

sample='sample_name'
outdir='output_dir'
#change work directory
cd $outdir/$sample/
#run the runWorkflow.py
$outdir/$sample/runWorkflow.py -m local -j 8


