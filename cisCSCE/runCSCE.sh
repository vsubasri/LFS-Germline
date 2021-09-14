#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.CSCE.err
#PBS -o log.CSCE.out
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -l walltime=248:00:00

module load R/3.4.4 bedtools/2.27.1 bedops/2.4.3 samtools/1.5

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/CSCE_meQTL.R --value Noob_beta_850k_ComBatProject_PEER_apcsgt
Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/CSCE_meQTL.R --value Noob_beta_450k_ComBatProject_PEER_apcsgt
Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/CSCE_meQTL.R --value Noob_beta_nci_ComBatProject_PEER_apcsgt

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/CSCE_overlap_discovery_validation.R
Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/CSCE_meQTL_EWAS.R


