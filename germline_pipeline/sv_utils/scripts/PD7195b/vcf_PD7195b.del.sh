#!/bin/sh
# Generated using Pype

#PBS -N vcf_PD7195b.del
#PBS -l nodes=1:ppn=1
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l walltime=240:00:00
#PBS -e /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/logs/PD7195b
#PBS -o /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/logs/PD7195b
#PBS -W depend=afterok:49992745:49992746

module load python/3.4.0
module load samtools/1.2

python /hpf/largeprojects/adam/local/sw/pype/2.3.1/vcf2tab/vcf2tab.py --vcf /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/normal_PD7195b.del.tab.vcf \
                                            --output /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/ \
                                            --bam /hpf/largeprojects/adam/projects/lfs/data/sanger/PD7195b.wgs/alignment/PD7195b.realigned-recalibrated.bam \
                                            --name normal_PD7195b.del.tab \
                                            --lowqual True

python /hpf/largeprojects/adam/local/sw/pype/2.3.1/vcf2tab/vcf2tab.py --vcf /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.del.tab.vcf \
                                            --output /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/ \
                                            --bam /hpf/largeprojects/adam/projects/lfs/lfs_germline/aligned_bams/NA12878/NA12878_20k.b37.bam \
                                            --name tumor_PD7195b.del.tab 