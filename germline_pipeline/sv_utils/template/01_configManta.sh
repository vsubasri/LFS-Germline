#!/bin/bash
module purge; module load python/2.7.10
#generate MantaWorkflow
nor='sample_name'
bam='sample_path'
outdir='output_dir'

/hpf/largeprojects/adam/local/sw/manta/0.20.2/bin/configManta.py \
--normalBam=${bam} \
--referenceFasta=/hpf/largeprojects/adam/local/reference/homosapiens/ucsc/hs37d5/fasta/hs37d5.fa \
--runDir=${outdir}/${nor}

