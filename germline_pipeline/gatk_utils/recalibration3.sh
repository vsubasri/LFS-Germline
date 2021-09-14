#!/bin/bash
#PBS -S /bin/bash
#PBS -N Recalibration
#PBS -e log/VQSR.v3.err
#PBS -o log/VQSR.v3.out
#PBS -l mem=48g
#PBS -l vmem=124g
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load bedtools/2.25.0 java/1.8.0_91

ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
interval_list=/hpf/largeprojects/adam/projects/lfs/resources/b37_wgs_calling_regions.interval_list
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-3.8/GenomeAnalysisTK.jar
db_snp=/hpf/largeprojects/adam/projects/lfs/resources/dbsnp_138.b37.vcf
hapmap=/hpf/largeprojects/adam/projects/lfs/resources/hapmap_3.3.b37.vcf
omni=/hpf/largeprojects/adam/projects/lfs/resources/1000G/1000G_omni2.5.b37.vcf
thouG=/hpf/largeprojects/adam/projects/lfs/resources/1000G/1000G_phase1.snps.high_confidence.b37.vcf
mills=/hpf/largeprojects/adam/projects/lfs/resources/Mills_and_1000G_gold_standard.indels.b37.vcf

TMPDIR=/localhd/$PBS_JOBID
BATCH=$(echo $(basename $INPUTFILE) | sed 's/genotyped_variants.\(.*\).g.vcf/\1/')
GATKDIR=$(dirname $INPUTFILE)

function snp_recalibration {
	java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T VariantRecalibrator \
		-R $ref_fasta \
		-input ${GATKDIR}genotyped_variants.v3.${BATCH}.vcf \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 $thouG \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $db_snp \
		-an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
		-mode SNP \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile ${GATKDIR}raw.SNPs.recal.v3.${BATCH}.vcf \
		-tranchesFile ${GATKDIR}raw.SNPs.v3.tranches.${BATCH} \
		-rscriptFile  ${GATKDIR}recal.plots.v3.${BATCH}.R \
		-nt 16

	java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T ApplyRecalibration \
		-R $ref_fasta \
		-input ${GATKDIR}genotyped_variants.v3.${BATCH}.vcf \
		-mode SNP \
		--ts_filter_level 99.5 \
		-recalFile ${GATKDIR}raw.SNPs.recal.v3.${BATCH}.vcf \
		-tranchesFile ${GATKDIR}raw.SNPs.v3.tranches.${BATCH} \
		-nt 40 \
		-o ${GATKDIR}SNPs.recalibrated.v3.${BATCH}.vcf
}

function indel_recalibration {
	java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T VariantRecalibrator \
		-R $ref_fasta \
		-input ${GATKDIR}SNPs.recalibrated.v3.${BATCH}.vcf \
		--maxGaussians 4 \
		-resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $db_snp \
		-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
		-mode INDEL \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile ${GATKDIR}raw.indels.recal.v3.${BATCH}.vcf \
		-tranchesFile ${GATKDIR}raw.indels.v3.tranches.${BATCH} \
		-rscriptFile ${GATKDIR}recal.plots.v3.${BATCH}.R \
		-nt 16TCH}

	java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk -T ApplyRecalibration \
		-R $ref_fasta \
		-input ${GATKDIR}SNPs.recalibrated.v3.${BATCH}.vcf \
		-mode INDEL \
		--ts_filter_level 99.0 \
		-recalFile ${GATKDIR}raw.indels.recal.v3.${BATCH}.vcf \
		-tranchesFile ${GATKDIR}raw.indels.v3.tranches.${BATCH} \
		-nt 40 \
		-o ${GATKDIR}indel.recalibrated.v3.${BATCH}.vcf
}

snp_recalibration
indel_recalibration
