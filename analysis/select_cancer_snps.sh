#!/bin/bash
#PBS -S /bin/bash
#PBS -N SelectVariants
#PBS -e select_variants.err
#PBS -o select_variants.out
#PBS -l mem=24g
#PBS -l vmem=24g
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 picard-tools/2.5.0 python/3.5.2

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=/localhd/$PBS_JOBID
input=/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/gatk/indel.recalibrated.vcf
output_vcf=/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/gatk/cancer_rsSNP.vcf
output_tab=/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/gatk/cancer_rsSNP.tab
fileKeep=/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/gatk/cancer_rsSNPids.txt

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx8G $gatk \
   -R ${ref_fasta} \
   -T SelectVariants \
   --variant ${input} \
   -o ${output_vcf} \
   -IDs ${fileKeep}

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx8G $gatk \
     -R ${ref_fasta} \
     -T VariantsToTable \
     -V ${output_vcf} \
     -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F ID -GF GT \
     -o ${output_tab}

python /hpf/largeprojects/adam/projects/lfs/lfs_germline/gatk/unfiltered_variants/complete/extract_probesXsnvs.py -s $output_tab

