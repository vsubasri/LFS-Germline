#!/bin/bash

#PBS -S /bin/bash
#PBS -N v3GenotypeGVCFs
#PBS -e log/genotypegvcfs.v3.err
#PBS -o log/genotypegvcfs.v3.out
#PBS -l mem=32g
#PBS -l vmem=124g
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91 

####################### GATK CombineGVCFs   ##################################################################################
###
###
# Take as input the per sample per chromosome gVCF files and produces the combined gvcf file which is
# meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs
# Reference: https://gatkforums.broadinstitute.org/gatk/discussion/10061/using-genomicsdbimport-to-consolidate-gvcfs-for-input-to-genotypegvcfs-in-gatk4
###
###

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=/localhd/$PBS_JOBID

BATCH=$(echo $(basename $INPUTLIST) | sed 's/input.batch\(.*\).list/\1/')
GATKDIR=$(dirname $INPUTLIST)

function CombineGVCFs {
    ##ImportGenomicsDB faster than CombineGVCFs
        java -jar -Djava.io.tmpdir=$TMPDIR -Xmx72G $gatk -T CombineGVCFs \
        --variant ${GATKDIR}input.${BATCH}.list \
        -R ${ref_fasta} \
        -o ${GATKDIR}combined_variants.${BATCH}.v3.g.vcf
}

####################### GATK GenotypeGVCFs  ##################################################################################
###
###
# This is the part that combines all the VCFs across samples to do the joint calling.
# This is a more practical aprroach of doing joint-calling than using the UnifiedGenotyper
# which relies on the BAM files.
# Reference for annotation options: https://software.broadinstitute.org/gatk/documentation/article?id=11069
###
###

function GenotypeGVCFs {
    java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T GenotypeGVCFs \
    	-R ${ref_fasta} \
    	-L ${interval_list} \
    	--interval_padding 100 \
    	-D ${db_snp} \
    	--variant ${GATKDIR}combined_variants.${BATCH}.v3.g.vcf \
    	-o ${GATKDIR}genotyped_variants.${BATCH}.v3.vcf \
    	--annotation InbreedingCoeff \
    	--annotation QualByDepth \
    	--annotation MappingQualityRankSumTest \
    	--annotation ReadPosRankSumTest \
    	--annotation FisherStrand 
}

CombineGVCFs
GenotypeGVCFs
