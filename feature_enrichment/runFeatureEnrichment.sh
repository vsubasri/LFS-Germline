#!/bin/bash
#PBS -S /bin/bash
#PBS -l vmem=4g
#PBS -l mem=4g
#PBS -N Feature_Enrichment
#PBS -e Feature_Enrichment.err
#PBS -o Feature_Enrichment.out

module load kentutils/302.1 bedtools/2.27.1

marks=(RNAseq H3K4me1 H3K4me3 H3K27me3 H3K9me3 H3K27ac H3K36me3 DNase H2A.Z H3K79me2)

wd=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment
enrichment=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/enrichment.sh
enrichment_mean=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/enrichment_mean.sh
infile=lfs_ageofonset_seed6543_S25_NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR.bed

if [ ! -d $wd/output ]; then
	mkdir $wd/output
fi

if [ ! -d $wd/output_mean ]; then
        mkdir $wd/output_mean
fi

cd $wd
for mark in ${marks[@]}
do
	if [ ! -d $mark ]; then
		mkdir $mark
	fi
	
	if [ ! -d $wd/output_mean/$mark ]; then
		mkdir $wd/output_mean/$mark
	fi

        if [ ! -d $wd/output/$mark ]; then
                mkdir $wd/output/$mark
        fi

	for i in E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046 E030;
	do 
		echo "$i $mark"

		## Take enrichment sum ##
#		if [[ ! -s $wd/output/$mark/predictive_$sample-$mark.bed ]] || [[ ! -s $wd/output/$mark/sampled_$sample-$mark.bed ]]; then
#			qsub -v sample=$i,mark=$mark,wd=$wd $enrichment	
#
		## Take enrichment mean ##
		
		if [[ ! -s $wd/output_mean/$mark/predictive_$sample-$mark.bed ]] || [[ ! -s $wd/output_mean/$mark/sampled_$sample-$mark.bed ]]; then
			qsub -v sample=$i,mark=$mark,wd=$wd,infile=$infile $enrichment_mean
		fi
	done
done
