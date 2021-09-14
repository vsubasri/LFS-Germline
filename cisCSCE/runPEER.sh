#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.PEER.err
#PBS -o log.PEER.out
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -l walltime=720:00:00

module load R/3.4.4

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

id=Noob_beta_nci_ComBatProject

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/PEER.R --value ${id}

