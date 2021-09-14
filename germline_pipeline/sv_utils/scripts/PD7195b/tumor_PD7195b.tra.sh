#!/bin/sh
# Generated using Pype

#PBS -N tumor_PD7195b.tra
#PBS -l nodes=1:ppn=1
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l walltime=240:00:00
#PBS -e /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/logs/PD7195b/ 
#PBS -o /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/logs/PD7195b/

module load python/3.4.0
module load samtools/1.2

export LC_ALL=C ;



/hpf/largeprojects/adam/local/sw/delly/0.7.1S/delly -t TRA \
                                                   -x human.hg19.excl.tsv \
                                                   -o /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra.tab.vcf \
                                                   -g /hpf/largeprojects/adam/local/reference/homosapiens/ucsc/hs37d5/fasta/hs37d5.fa /hpf/largeprojects/adam/projects/lfs/lfs_germline/aligned_bams/NA12878/NA12878_20k.b37.bam \
                                                   -p /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra.tab.pe.metrics.txt
#                                                   -r /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/tumor_PD7195b.tra.tab.sr.metrics.txt

python /hpf/largeprojects/adam/local/sw/pype/2.3.1/bam2vcf/generateBAM.py --metrics /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra.tab.pe.metrics.txt \
                                                --bamfile /hpf/largeprojects/adam/projects/lfs/lfs_germline/aligned_bams/NA12878/NA12878_20k.b37.bam \
                                                --output  /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_PE.bam
samtools sort /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_PE.bam /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_temp
mv /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_temp.bam /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_PE.bam
samtools index /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra_PE.bam
#python /hpf/largeprojects/adam/local/sw/pype/2.3.1/bam2tab/generateBAM.py --metrics /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/tumor_PD7195b.tra.tab.sr.metrics.txt \
#                                                --bamfile /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/tumor_PD7195b.tra.tab \
#                                                --output  /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/tumor_PD7195b.tra_SR.bam
rm /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/data/PD7195b/tumor_PD7195b.tra.tab.pe.metrics.txt
# rm tumor_PD7195b.tra.tab.sr.metrics.txt

