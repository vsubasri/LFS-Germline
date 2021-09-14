#!/bin/bash
#PBS -S /bin/bash
#PBS -e Preprocessing.err
#PBS -o Preprocessing.out
#PBS -l vmem=72g
#PBS -l mem=72g
#PBS -l walltime=448:00:00

module load R/3.4.4

outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/
id=Noob_beta_nci

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/preprocess_beta.R

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/batch_correction_ComBat.R --value ${id}_ComBatProject --infile ${id}.rds --outdir ${outdir}

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/PEER.R --value ${id}_ComBatProject

